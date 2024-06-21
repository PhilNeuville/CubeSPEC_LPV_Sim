from functions import *
import matplotlib.pyplot as plt
import inspect
import imageio.v2 as imageio

# Pulsation time serie
rv, pulsations = Get_pulsations(pulsation_dir)                      # extract pulsations info
rev_pul = Line_invert(pulsations)                                   # reverse pulsations
mean_profile = Mean_pul(rev_pul)                                    # compute mean profile
normed_pul = Line_invert(Norm_pul(rv, rev_pul, mean_profile))       # normalise pulsations

for i in range(np.shape(normed_pul)[1]):
    plt.plot(rv, pulsations[:,i])
    # plt.plot(rv, pulsations[:, i])
plt.show()

wvl = (rv / (cc * 1e-3) + 1) * 5200

broadPul = np.zeros_like(pulsations)
for pul in range(np.shape(pulsations)[1]):
    broadPul[:, pul] = Apply_Rotation(wvl, pulsations[:, pul], 100)

for j in range(np.shape(broadPul)[1]):
    plt.plot(rv, normed_pul[:, j])
    plt.plot(rv, broadPul[:, j])
plt.show()

test_rotBroad = Apply_Rotation(wvl, pulsations[:, 0], 100)

exit()

# Null profile
# _, null_profile = Get_null_profile(null_profile_dir)
# rev_null = Line_invert(null_profile)
# normed_null = rev_null / np.trapz(rev_null, rv)

# zz = np.zeros(100)
# new = np.hstack((zz, rev_pul[:, 0], zz))
# quotient, remainder = signal.deconvolve(new, mean_profile)    # deconvolve thingy
# sig = signal.convolve(mean_profile, quotient, 'same') + remainder[100:-100]

# quotient, remainder = signal.deconvolve(normed_pul[:, 0], normed_null) #---> still not working

# # Plotting tests
# fig, ax = plt.subplots(figsize = (6,4), dpi = 150)
#
# for i in range(np.shape(pulsations)[1]):
#
#     # ax.plot(rv, pulsations[:, i], alpha = 0.2)
#     ax.plot(rv, rev_pul[:, i], alpha = 0.1)
#     # ax.plot(rv, normed_pul[:, i], alpha = 0.1)
#
# # ax.plot(wvl, normed_pul[:, 0])
#
# ax.plot(rv, mean_profile, label ='inverted mean profile')
# ax.plot(rv, rev_null, label = 'null profile')
# # ax.plot(rv, normed_null, label = 'normalised null profile')
# # ax.plot(rv, sig, color = 'red', label='recovered signal')
# # ax.plot(rv, remainder[100:-100], label ='remainder')
# # ax.plot(rv, rev_pul[:, 0], color = 'green', ls = '--')
# # ax.plot(rv, quotient, color = 'green', ls = '--')
#
# plt.legend(loc='best')
# plt.tight_layout()
# plt.show()

# Test to save multiple files in a folder
#give_outputs(output_path, 'Test2', rv, normed_pul)

#print(inspect.getsource(pyasl.instrBroadGaussFast))


# print(np.shape(pulsations))
#
# test = np.zeros((2, 3))
# print(test)
# print(np.shape(test))