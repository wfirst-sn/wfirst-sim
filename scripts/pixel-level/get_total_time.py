from numpy import *
import commands
import matplotlib.pyplot as plt
from pixel_level_ETC import *
from copy import deepcopy


slew_time = 42.
generate_data = 0
SN_counts = [69, 208, 403, 223, 327, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136]
#SN_counts = [69] + [400]*6 + [80]*9

runs = [#(0.15, 1, "test", 4.0, 0),
        #(0.15, 1, "AiryObstruct", 4.0, 0),
        #(0.15, 1, "WebbPSF", 4.0, 0),
        #(0.075, 2, "WebbPSF", 4.0, 0),
        #(0.075, 2, "WebbPSF", 3, 0),
        #(0.15, 1, "WebbPSF", 4, 6e-19),
        (0.075, 2, "WebbPSF", 4, 6e-19),
        
        #(0.15, 1, "WebbPSF", 3.),
        #(0.075, 2, "WebbPSF", 3.)
        #(0.15, 1, "NIC"),
        #(0.075, 2, "NIC")]:
    ]

offset_i = 0
offset_j = 10.
dark_current = 0.01
read_noise = None

if 0:
    #runs = [(0.15, 1, "WebbPSF", read_noise_floor, 0) for read_noise_floor in arange(0.0, 5.01, 0.5)]
    runs = [(0.075, 1, "WebbPSF", read_noise_floor, 0) for read_noise_floor in arange(0.0, 5.01, 0.5)]

for pixel_scale, slice_in_pixels, PSF_source, read_noise_floor, gal_flamb in runs:

    PSFs = initialize_PSFs(scales = [int(round(pixel_scale/0.005))], PSF_source = PSF_source)
    redshifts = arange(0.15, 1.66, 0.1)
    args = {"gal_flamb": lambda x:gal_flamb, "pixel_scale": pixel_scale, "slice_in_pixels": slice_in_pixels, "dark_current": dark_current, "mdl": 'hsiao', "PSFs": PSFs, "IFURfl": "IFU_R_Content.txt", "min_wave": 4000., "read_noise": read_noise, "read_noise_floor": read_noise_floor, "offset_i": offset_i, "offset_j": offset_j}

    if 0:
        args["key1"] = "obs_frame"
        args["key2"] = (10200, 12850)
        StoNs = [3.5, 6.0, 10]
    elif 1:
        args["key1"] = "rest_frame_mean_S/N"
        args["key2"] = (5000, 6000)
        StoNs = [3.5, 6.0, 10]
    elif 0:
        args["key1"] = "rest_frame_band_S/N"
        args["key2"] = (5000, 6000)
        StoNs = [14.5, 25, 41.5]
    
    total_time = 0.
    j = 0
    plt.figure(figsize = (60,40))

    for SN_count, redshift in zip(SN_counts, redshifts):
        # eight at S/N 3.5 + three at S/N 6 + one at S/N 10 = 8 + 3*2 + 4 = 18
        exptime35 = solve_for_exptime(S_to_N = StoNs[0], redshift = redshift, **args)
        exptime6 = solve_for_exptime(S_to_N = StoNs[1], redshift = redshift, **args)
        exptime10 = solve_for_exptime(S_to_N = StoNs[2], redshift = redshift, **args)
        
        time_at_this_redshift_per_SN = ((exptime35 + slew_time)*8 + (exptime6 + slew_time)*3 + (exptime10 + slew_time)*1)
        total_time += time_at_this_redshift_per_SN*SN_count


        toprint = [redshift, "%.1f: %.1f" % (StoNs[0], exptime35), "%.1f: %.1f" % (StoNs[1], exptime6), "%.1f: %.1f" % (StoNs[2], exptime10), SN_count, "%.1f" % time_at_this_redshift_per_SN,
                   "Days: %.1f" % (time_at_this_redshift_per_SN*SN_count/86400.), "%.1f" % (total_time/86400.), "offset %i %i" % (args["offset_i"], args["offset_j"]),
                   "PSF: " + PSF_source, str(pixel_scale) + '"', str(slice_in_pixels) + " pixel/slice", "RN floor: " + str(read_noise_floor), "gal_flamb: %.1g" % gal_flamb, args["key1"] + " " + str(args["key2"])]

        toprint = [str(item) for item in toprint]
        print '\t'.join(toprint)


        phase_times = [(-6, exptime35),
                       (-2, exptime6),
                       (2, exptime10),
                       (5, exptime35),
                       (8, exptime35),
                       (11, exptime35),
                       (14, exptime35),
                       (17, exptime35),
                       (20, exptime35)]


        if generate_data:

            for i, item in enumerate(phase_times):
                plt.subplot(len(SN_counts), len(phase_times), i+1 + j*len(phase_times))
                new_args = deepcopy(args)
                del new_args["key1"]
                del new_args["key2"]
                ETC_dict = get_spec_with_err(redshift = redshift, exp_time = item[1], phase = item[0], show_plots = 0, **new_args)
                
                
                plt.plot(ETC_dict["obs_waves"], ETC_dict["spec_S/N"], '.', color = 'k')
                #plt.plot(ETC_dict["obs_waves"], ETC_dict["f_lamb_SN"], color = 'k')
                #ylim = plt.ylim()
                #plt.fill_between(ETC_dict["obs_waves"], ETC_dict["f_lamb_SN"]*(1 - 1./ETC_dict["spec_S/N"]), ETC_dict["f_lamb_SN"]*(1 + 1./ETC_dict["spec_S/N"]), color = (0.8, 0.8, 0.8))
                #plt.ylim(ylim)

        j += 1
    
    if generate_data:
        plt.savefig("output/survey_grid_" + ETC_dict["plt_root"] + ".pdf", bbox_inches = 'tight')
        plt.close()
    print
    print

    
