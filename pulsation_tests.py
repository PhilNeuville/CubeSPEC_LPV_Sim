from functions import *
import matplotlib.pyplot as plt
import inspect

# Pulsation time serie
wvl, pulsations = get_pulsations(pulsation_dir)        # extract pulsations info
rev_pul = line_invert(pulsations)                        # reverse pulsations
mean_profile = mean_pul(rev_pul)                         # compute mean profile
normed_pul = norm_pul(wvl, rev_pul, mean_profile)        # normalise pulsations

# Null profile
_, null_profile = get_null_profile(null_profile_dir)
rev_null = line_invert(null_profile)
normed_null = rev_null / np.trapz(rev_null, wvl)

# zz = np.zeros(100)
# new = np.hstack((zz, rev_pul[:, 0], zz))
# quotient, remainder = signal.deconvolve(new, mean_profile)    # deconvolve thingy
# sig = signal.convolve(mean_profile, quotient, 'same') + remainder[100:-100]

# quotient, remainder = signal.deconvolve(normed_pul[:, 0], normed_null) #---> still not working

if __name__ == '__main__':

    # Plotting tests
    fig, ax = plt.subplots(figsize = (6,4), dpi = 150)

    for i in range(np.shape(pulsations)[1]):

        ax.plot(wvl, pulsations[:, i], alpha = 0.2)
        # ax.plot(wvl, rev_pul[:, i], alpha = 0.1)
        # ax.plot(wvl, normed_pul[:, i], alpha = 0.1)

    # ax.plot(wvl, normed_pul[:, 0])

    # ax.plot(wvl, mean_profile, label ='inverted mean profile')
    # ax.plot(wvl, rev_null, label = 'null profile')
    # ax.plot(wvl, normed_null, label = 'normalised null profile')
    # ax.plot(wvl, sig, color = 'red', label='recovered signal')
    # ax.plot(wvl, remainder[100:-100], label ='remainder')
    # ax.plot(wvl, rev_pul[:, 0], color = 'green', ls = '--')
    # ax.plot(wvl, quotient, color = 'green', ls = '--')

    #plt.legend(loc='best')
    plt.tight_layout()
    plt.show()

# Test to save multiple files in a folder
#give_outputs(output_path, 'Test2', wvl, normed_pul)

#print(inspect.getsource(pyasl.instrBroadGaussFast))
