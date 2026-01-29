import torch
from torch.nn.functional import max_pool3d, conv3d, group_norm, leaky_relu, interpolate, relu, batch_norm
import numpy as np
from SuperSynth.utils import MRIread

##############

class frugal_models():
    def __init__(self, model_file, MNIqcseg_file, device):

        print('Preparing model and loading weights')
        self.device = device
        if device == 'cpu':
            cp = torch.load(model_file, map_location=torch.device('cpu'))
        else:
            cp = torch.load(model_file)

        # SuperSynth part
        self.ssynth_bbone = cp['backbone_state_dict']
        for key in self.ssynth_bbone.keys():
            self.ssynth_bbone[key] = self.ssynth_bbone[key].to(device)
        self.ssynth_final_conv_names = ['reg', 'seg', 'T1', 'T2', 'FLAIR', 'LP', 'LW', 'RP', 'RW']
        self.ssynth_final_conv_weight = {}
        self.ssynth_final_conv_bias = {}
        for name in self.ssynth_final_conv_names:
            self.ssynth_final_conv_weight[name] = cp[name + '_state_dict']['weight'].to(device)
            self.ssynth_final_conv_bias[name] = cp[name + '_state_dict']['bias'].to(device)

        # AutoQC part (Billot et al, PNAS, 2023)
        for key in cp.keys():
            if key.endswith('_state_dict') is False:
                setattr(self, key, cp[key].to(device))
        labels_segmentation = np.array(
            [0, 2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 24, 26, 41, 42, 43, 44, 46, 47, 49, 50, 51, 52,
             53, 54, 58, 77, 85, 99, 819, 820, 821, 822, 843, 844, 865, 866, 869, 870, 901, 902, 906, 907, 908, 909,
             911, 912, 914, 915, 916, 930], dtype=np.int32)
        labels_qc = np.array(
            [0, 1, 2, 3, 3, 4, 4, 6, 1, 7, 7, 3, 3, 5, 8, 8, 0, 1, 1, 2, 3, 3, 4, 4, 6, 1, 7, 7, 8, 8, 1, 1, 1, 0, 1, 1,
             1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=np.int32)
        self.lut = torch.zeros(1000, dtype=torch.int32, device=device)
        for i in range(len(labels_qc)):
            self.lut[labels_segmentation[i]] = labels_qc[i]
        self.onehotlut = torch.eye(9, dtype=torch.float32, device=device)

        # Clean up
        del cp
        torch.cuda.empty_cache()

        # Also, prepare stuff we need to inpaint missing regions for QC module
        MNIqcseg, self.aff_mni = MRIread(MNIqcseg_file)
        self.MNIqcseg = torch.tensor(MNIqcseg, device=device, dtype=torch.int32)
        ii = torch.arange(MNIqcseg.shape[0], device=device, dtype=torch.float32,)
        jj = torch.arange(MNIqcseg.shape[1], device=device, dtype=torch.float32,)
        kk = torch.arange(MNIqcseg.shape[2], device=device, dtype=torch.float32,)
        self.ii, self.jj, self.kk = torch.meshgrid(ii, jj, kk, indexing='ij')


    def ssynth_inference(self, input):

        enc_feat_maps = []
        for l in range(5):
            x = input if l == 0 else max_pool3d(enc_feat_maps[-1], 2)
            x = conv3d(x, self.ssynth_bbone['encoders.' + str(l) + '.basic_module.conv1.weight'],
                       self.ssynth_bbone['encoders.' + str(l) + '.basic_module.conv1.bias'])
            x = conv3d(x, self.ssynth_bbone['encoders.' + str(l) + '.basic_module.conv2.conv.weight'], padding=1)
            x = group_norm(x, 8, weight=self.ssynth_bbone['encoders.' + str(l) + '.basic_module.conv2.groupnorm.weight'],
                           bias=self.ssynth_bbone['encoders.' + str(l) + '.basic_module.conv2.groupnorm.bias'])  # , eps=1e-05)
            x = leaky_relu(x, negative_slope=0.01, inplace=True)
            x = conv3d(x, self.ssynth_bbone['encoders.' + str(l) + '.basic_module.conv3.conv.weight'], padding=1)
            x = group_norm(x, 8, weight=self.ssynth_bbone['encoders.' + str(l) + '.basic_module.conv3.groupnorm.weight'],
                           bias=self.ssynth_bbone['encoders.' + str(l) + '.basic_module.conv3.groupnorm.bias'])  # , eps=1e-05)
            # next line is a bit inefficient but saves memory!
            x += conv3d(input if l == 0 else max_pool3d(enc_feat_maps[-1], 2),
                        self.ssynth_bbone['encoders.' + str(l) + '.basic_module.conv1.weight'],
                        self.ssynth_bbone['encoders.' + str(l) + '.basic_module.conv1.bias'])
            x = leaky_relu(x, negative_slope=0.1, inplace=True)
            enc_feat_maps.append(x)
            del x
            torch.cuda.empty_cache()

        for l in range(4):
            # Here's the real memory bottleneck that we overcome with ugly code
            idx = enc_feat_maps[-2].shape[1]
            nc = 24  # how many channels at a time
            enc_feat_maps[-2] = conv3d(enc_feat_maps[-2],
                                       self.ssynth_bbone['decoders.' + str(l) + '.basic_module.SingleConv1.conv.weight'][:, :idx],
                                       padding=1)
            newshape = enc_feat_maps[-2].shape[2:]
            for c in range(0, enc_feat_maps[-1].shape[1], nc):
                enc_feat_maps[-2] += conv3d(interpolate(enc_feat_maps[-1][:, c:c + nc], size=newshape, mode='nearest'),
                                            self.ssynth_bbone['decoders.' + str(l) + '.basic_module.SingleConv1.conv.weight'][:, idx + c:idx + nc + c],
                                            padding=1)

            del enc_feat_maps[-1]  # now, enc_feat_maps[-2] is enc_feat_maps[-1]
            torch.cuda.empty_cache()
            enc_feat_maps[-1] = group_norm(enc_feat_maps[-1], 8, weight=self.ssynth_bbone[
                'decoders.' + str(l) + '.basic_module.SingleConv1.groupnorm.weight'],
                                           bias=self.ssynth_bbone['decoders.' + str(l) + '.basic_module.SingleConv1.groupnorm.bias'],
                                           eps=1e-05)
            enc_feat_maps[-1] = leaky_relu(enc_feat_maps[-1], negative_slope=0.01, inplace=True)
            enc_feat_maps[-1] = conv3d(enc_feat_maps[-1],
                                       self.ssynth_bbone['decoders.' + str(l) + '.basic_module.SingleConv2.conv.weight'],
                                       padding=1)
            enc_feat_maps[-1] = group_norm(enc_feat_maps[-1], 8, weight=self.ssynth_bbone[
                'decoders.' + str(l) + '.basic_module.SingleConv2.groupnorm.weight'],
                                           bias=self.ssynth_bbone['decoders.' + str(l) + '.basic_module.SingleConv2.groupnorm.bias'],
                                           eps=1e-05)
            enc_feat_maps[-1] = leaky_relu(enc_feat_maps[-1], negative_slope=0.01, inplace=True)

        outputs = {}
        for name in self.ssynth_final_conv_names:
            outputs[name] = conv3d(enc_feat_maps[0], self.ssynth_final_conv_weight[name], self.ssynth_final_conv_bias[name])

        del enc_feat_maps
        torch.cuda.empty_cache()

        return outputs


    def inpaint_seg_as_needed_for_qc(self, seg, mode, seg_aff, mni_affine):

        if (mode == 'invivo') or (mode == 'exvivo'):  # nothing to do!
            seg_inpainted = seg
        else:
            A = np.linalg.inv(seg_aff) @ np.linalg.inv(mni_affine.detach().cpu().numpy().astype(np.float32)) @ self.aff_mni
            ii2 = (A[0, 0] * self.ii + A[0, 1] * self.jj + A[0, 2] * self.kk + A[0, 3]).round().to(seg.dtype)
            jj2 = (A[1, 0] * self.ii + A[1, 1] * self.jj + A[1, 2] * self.kk + A[1, 3]).round().to(seg.dtype)
            kk2 = (A[2, 0] * self.ii + A[2, 1] * self.jj + A[2, 2] * self.kk + A[2, 3]).round().to(seg.dtype)
            ok = (ii2 >= 0) & (jj2 >= 0) & (kk2 >= 0) & (ii2 <= (seg.shape[0] - 1)) & (jj2 <= (seg.shape[1] - 1)) & (kk2 <= (seg.shape[2] - 1))
            vals = torch.zeros_like(ii2)
            vals[ok] = seg[ii2[ok], jj2[ok], kk2[ok]]
            seg_warped = vals.reshape(self.MNIqcseg.shape)
            if (mode == 'left-hemi') or (mode == 'right-hemi'):
                seg_warped = torch.maximum(seg_warped, seg_warped.flip([0]))
            seg_inpainted = torch.zeros_like(self.MNIqcseg)
            mask = (self.MNIqcseg == 3) | (self.MNIqcseg == 4) | (self.MNIqcseg == 5)
            seg_inpainted[mask] = self.MNIqcseg[mask]
            mask = seg_warped > 0
            seg_inpainted[mask] = seg_warped[mask].to(seg_inpainted.dtype)


        return seg_inpainted

    def qc_inference(self, input, mode='invivo', input_aff=None, mni_affine=None):

        # Map to QC labels and inpaint if needed
        seg = self.inpaint_seg_as_needed_for_qc(self.lut[input], mode, input_aff, mni_affine)

        # Pad to 224x224x224 (ask Benjamin...)
        w = torch.where(seg > 0)
        seg = seg[w[0].min():(1 + w[0].max()), w[1].min():(1 + w[1].max()), w[2].min():(1 + w[2].max())]
        seg_padded = torch.zeros([224, 224, 224], device=self.device, dtype=torch.int32)
        in_shape = np.array(seg.shape)
        in_start = np.maximum((in_shape - 224) // 2, 0)
        out_start = np.maximum((224 - in_shape) // 2, 0)
        copy_size = np.minimum(in_shape, 224)
        in_slices = tuple(slice(in_start[i], in_start[i] + copy_size[i]) for i in range(3))
        out_slices = tuple(slice(out_start[i], out_start[i] + copy_size[i]) for i in range(3))
        seg_padded[out_slices] = seg[in_slices]
        onehot = self.onehotlut[seg_padded]

        # Level 1
        inputLevel = onehot.permute([3,0,1,2])[None, ...]
        last_tensor = conv3d(inputLevel, self.conv_downarm_0_0_kernel, bias=self.conv_downarm_0_0_bias, padding='same')
        last_tensor = relu(last_tensor, inplace=True)
        last_tensor = conv3d(last_tensor, self.conv_downarm_0_1_kernel, bias=self.conv_downarm_0_1_bias, padding='same')
        last_tensor += relu(conv3d(inputLevel, self.conv_expand_0_kernel, bias=self.conv_expand_0_bias, padding='same'), inplace=True)
        last_tensor = relu(last_tensor, inplace=True)
        last_tensor = batch_norm(last_tensor, running_mean=self.bn_0_moving_mean, running_var=self.bn_0_moving_variance,
                         weight=self.bn_0_gamma, bias=self.bn_0_beta, training=False, momentum=0, eps=0.001)
        last_tensor = max_pool3d(last_tensor, 2)

        # Level 2
        inputLevel = last_tensor
        last_tensor = conv3d(last_tensor, self.conv_downarm_1_0_kernel, bias=self.conv_downarm_1_0_bias, padding='same')
        last_tensor = relu(last_tensor, inplace=True)
        last_tensor = conv3d(last_tensor, self.conv_downarm_1_1_kernel, bias=self.conv_downarm_1_1_bias, padding='same')
        last_tensor += relu(conv3d(inputLevel, self.conv_expand_1_kernel, bias=self.conv_expand_1_bias, padding='same'), inplace=True)
        last_tensor = relu(last_tensor, inplace=True)
        last_tensor = batch_norm(last_tensor, running_mean=self.bn_1_moving_mean, running_var=self.bn_1_moving_variance,
                         weight=self.bn_1_gamma, bias=self.bn_1_beta, training=False, momentum=0, eps=0.001)
        last_tensor = max_pool3d(last_tensor, 2)

        # Level 3
        inputLevel = last_tensor
        last_tensor = conv3d(last_tensor, self.conv_downarm_2_0_kernel, bias=self.conv_downarm_2_0_bias, padding='same')
        last_tensor = relu(last_tensor, inplace=True)
        last_tensor = conv3d(last_tensor, self.conv_downarm_2_1_kernel, bias=self.conv_downarm_2_1_bias, padding='same')
        last_tensor += relu(conv3d(inputLevel, self.conv_expand_2_kernel, bias=self.conv_expand_2_bias, padding='same'), inplace=True)
        last_tensor = relu(last_tensor, inplace=True)
        last_tensor = batch_norm(last_tensor, running_mean=self.bn_2_moving_mean, running_var=self.bn_2_moving_variance,
                         weight=self.bn_2_gamma, bias=self.bn_2_beta, training=False, momentum=0, eps=0.001)
        last_tensor = max_pool3d(last_tensor, 2)

        # Level 4
        inputLevel = last_tensor
        last_tensor = conv3d(last_tensor, self.conv_downarm_3_0_kernel, bias=self.conv_downarm_3_0_bias, padding='same')
        last_tensor = relu(last_tensor, inplace=True)
        last_tensor = conv3d(last_tensor, self.conv_downarm_3_1_kernel, bias=self.conv_downarm_3_1_bias, padding='same')
        last_tensor += relu(conv3d(inputLevel, self.conv_expand_3_kernel, bias=self.conv_expand_3_bias, padding='same'), inplace=True)
        last_tensor = relu(last_tensor, inplace=True)
        last_tensor = batch_norm(last_tensor, running_mean=self.bn_3_moving_mean, running_var=self.bn_3_moving_variance,
                         weight=self.bn_3_gamma, bias=self.bn_3_beta, training=False, momentum=0, eps=0.001)
        last_tensor = max_pool3d(last_tensor, 2)

        # Final convolutions and output
        last_tensor = conv3d(last_tensor, self.final_conv_0_kernel, bias=self.final_conv_0_bias, padding='same')
        last_tensor = relu(last_tensor, inplace=True)
        last_tensor = conv3d(last_tensor, self.final_conv_1_kernel, bias=self.final_conv_1_bias, padding='same')
        last_tensor = relu(last_tensor, inplace=True)
        qc_scores = last_tensor.mean(dim=[0,2,3,4])[1:].clip(0 ,1)

        # Mask cerebellum and brainstem scores if needed
        if (mode == 'left-hemi') or (mode == 'right-hemi') or (mode == 'cerebrum'):
            qc_scores[3] = qc_scores[4] = 0

        return qc_scores



