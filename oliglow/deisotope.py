import ms_deisotope as ms_ditp
import polars as pl
import importlib
import utils
import averageine

# get_mono from clr_loader defines where to get mono from
from clr_loader import get_mono

# Use the mono from the mono installation in the homebrew directory and not the default (the default is only for the intel acrchitecture)
rt = get_mono(libmono="/opt/homebrew/Cellar/mono/6.12.0.206/lib/libmono-2.0.1.dylib")

# import clr
import pythonnet

pythonnet.get_runtime_info()

MASS_PROTON = 1.00727646677  # Replace this by the charge of the mass carrier if using other than protonated ions!


def get_additional_bunch_info(bunch):
    # Get the relevant precursor information and MS2 scan number from the bunch object
    time = bunch.scan_time  # (in minutes)
    # ms2_scan_number = bunch.index + 1
    ms2_scan_number = int(bunch.scan_id.split("scan=")[-1])

    precursor_mz = bunch.precursor_information.mz
    ms1_precursor_scan_number = int(
        bunch.precursor_information.precursor_scan_id.split("scan=")[-1]
    )  # TODO: This does not match what freestyle does, check this!

    if type(bunch.precursor_information.charge) == type(1):
        precursor_charge = bunch.precursor_information.charge * bunch.polarity
        precursor_neutral_mass = (
            abs(precursor_mz * precursor_charge) + abs(precursor_charge) * MASS_PROTON
        )  # This reporrts a value close to the most_abundant_mass since
    else:
        precursor_charge = 0  # Sometimes the precursor charge is not known by the instrument, in that case, set it to 0!
        precursor_neutral_mass = 0.0

    ms1_intensity = bunch.precursor_information.intensity

    # Gives the osilation window target, lower and upper bounds of the MS1 preursor scan!
    isolation_window_target = bunch.isolation_window.target
    isolation_window_lower_bound = bunch.isolation_window.lower_bound
    isolation_window_upper_bound = bunch.isolation_window.upper_bound

    return pl.DataFrame(
        {
            "time": time,
            "ms2_scan_number": ms2_scan_number,
            "precursor_mz": precursor_mz,
            "precursor_neutral_mass": precursor_neutral_mass,
            "ms1_intensity": ms1_intensity,
            "ms1_precursor_scan_number": ms1_precursor_scan_number,
            "precursor_charge": precursor_charge,
            "isolation_window_target": isolation_window_target,
            "isolation_window_lower": isolation_window_lower_bound,
            "isolation_window_upper": isolation_window_upper_bound,
        }
    )


def re_deisotope_get_additional_bunch_info(bunch,precursor = None,re_deisotope_success = False):
    # Get the relevant precursor information and MS2 scan number from the bunch object

    if precursor is None:
        precursor = bunch.precursor_information

    time = bunch.scan_time  # (in minutes)
    # ms2_scan_number = bunch.index + 1
    ms2_scan_number = int(bunch.scan_id.split("scan=")[-1])

    ms1_precursor_scan_number = int(
        bunch.source.find_previous_ms1(bunch.index).scan_id.split("scan=")[-1]
        #bunch.precursor_information.precursor_scan_id.split("scan=")[-1]
    )  # TODO: This does not match what freestyle does, check this!

    # if re_deisotope_success:
    precursor_mz = precursor.extracted_mz
    ms1_intensity = precursor.extracted_intensity
    precursor_neutral_mass = precursor.extracted_neutral_mass

    if type(precursor.extracted_charge) == type(1) or type(precursor.extracted_charge) == type(-1):
        precursor_charge = precursor.extracted_charge
        # precursor_neutral_mass = precursor.extracted_neutral_mass
    else:
        if type(precursor.charge) == type(1) or type(precursor.charge) == type(-1):
            precursor_charge = precursor.charge  # Sometimes the precursor charge is not known by the instrument, in that case, set it to 0!
        else:
            precursor_charge = 0
        # precursor_neutral_mass = precursor.neutral_mass
    # else:
    #     precursor_mz = precursor.mz
    #     ms1_intensity = precursor.intensity
        
    #     if type(precursor.charge) == type(1) or type(precursor.charge) == type(-1):
    #         precursor_charge = precursor.charge
    #         precursor_neutral_mass = precursor.neutral_mass
    #     else:
    #         precursor_charge = 0  # Sometimes the precursor charge is not known by the instrument, in that case, set it to 0!
    #         precursor_neutral_mass = precursor.neutral_mass

    # Gives the osilation window target, lower and upper bounds of the MS1 preursor scan!
    isolation_window_target = bunch.isolation_window.target
    isolation_window_lower_bound = bunch.isolation_window.lower_bound
    isolation_window_upper_bound = bunch.isolation_window.upper_bound

    # TODO: Also read the UV intensity of the scan!

    return pl.DataFrame(
        {
            "time": time,
            "ms2_scan_number": ms2_scan_number,
            "precursor_mz": precursor_mz,
            "precursor_neutral_mass": precursor_neutral_mass,
            "ms1_intensity": float(ms1_intensity),
            "ms1_precursor_scan_number": ms1_precursor_scan_number,
            "ms1_deisotope_success": re_deisotope_success,
            "precursor_charge": precursor_charge,
            "isolation_window_target": isolation_window_target,
            "isolation_window_lower": isolation_window_lower_bound,
            "isolation_window_upper": isolation_window_upper_bound,
        }
    )


def deconvolute_ms2_bunch(
    bunch,
    parameters=dict(),
    re_deistope_precursor=True,
    trust_precursor_scan_number=False,
    re_deisotoping_parameters=dict(),
    return_data_frame=False
):
    # Deisotope a specific given bunch object!

    # Parameters

    if re_deisotoping_parameters == {}:
        re_deisotoping_parameters = parameters

    # TODO: Properly pass the re_deisotoping_parameters to the function bunch.precursor_information.find_monoisotopic_peak

    # Use the RNA averagine WITH the backbone!
    if "avgine" in parameters:
        avgine = parameters["avgine"]
    else:
        avgine = averageine.averagine_rna_with_backbone

    # The minimum score between the theo and exp spectra using MSDeconVFitter which is accepted!
    # TODO: Take care of this score, there should be way to estimate this score!#
    # Or use multiple passes, if the score is too large, i.e. the number of peaks after desiotoping are less than \alpha (num peaks), increase this value!
    if "min_score" in parameters:
        min_score = parameters["min_score"]
    else:
        min_score = 150.0

    if "mass_error_tol" in parameters:
        mass_error_tol = parameters["mass_error_tol"]
    else:
        mass_error_tol = 0.02
    # The Absolute error tolerance between the theoretical m/z and the exp m/z which is accepted!
    # perhaps this should be < 1./max(charge) for a reasonable value.

    scorer = ms_ditp.MSDeconVFitter(min_score, mass_error_tol)
    # The parameters are intensity scale depenedent.
    # We will need to find this dependant on the intensity of the scans!

    if "truncate_after" in parameters:
        truncate_after = parameters["truncate_after"]
    else:
        truncate_after = 0.8  # For MS1 and 0.8 for MS2 as per the recommendation of the ms_ditp author!
    # See the discussion in ms_isotope docs
    # SENSITIVE TO THIS PARAMETER! EXTREMELY SENSETIVE! Reduce this if high mass ions are missing!
    # The authors propose to use `incremental_truncation` if the latter is the case! Explore this later!
    # Play with the parameter Truncate_after! Perhaps also reduce the minimum score of the scorer, such that multiple peaks with problems are coalesced into one!

    if "scale" in parameters:
        scale = parameters["scale"]
    else:
        scale = "sum"  # What to do with the intensity of the different peaks!

    # If this paramter is increased then more deconvoluted peaks are selected, keep this small preferrentially (0 or 1)
    if "max_missed_peaks" in parameters:
        max_missed_peaks = parameters["max_missed_peaks"]
    else:
        max_missed_peaks = 0

    if "error_tol" in parameters:
        error_tol = parameters["error_tol"]
    else:
        #error_tol = 2e-5  # Also the default value! PPM error tolerance to use to match experimental to theoretical peaks, defalut = 2e-5
        error_tol = 2e-6

    if "scale_method" in parameters:
        scale_method = parameters["scale_method"]
    else:
        # scale_method = "max" #How to scale intensities when comparing theoretical to exp isotopic distributions, default = "sum"
        scale_method = "sum"  # How to scale intensities when comparing theoretical to exp isotopic distributions, default = "sum"

    # The maximum considered charge cannot be greater than that of the MS1 precursor charge!
    # Be careful of the polarity of the charge! If the charges are negative, the charge range needs to be supplied with negative sign!
    if "charge_range" in parameters:
        charge_range = parameters["charge_range"]
    else:
        if type(bunch.precursor_information.charge) == type(1):
            charge_range = (
                bunch.polarity,
                bunch.precursor_information.charge * bunch.polarity,
            )
        else:
            charge_range = (bunch.polarity, bunch.polarity * 30) #TODO: Estimate this from the sequence length
            #IMP: This number should be large enough to cover the charge states of the precursors!

    if "minimum_intensity" in parameters:
        minimum_intensity = parameters["minimum_intensity"]
    else:
        # minimum_intensity = 5. #Defualt = 5, ignore peaks below this intensity! Modify this based on the speatra!
        minimum_intensity = min(
            [bunch.peak_set[i].intensity for i in range(len(bunch.peak_set))]
        )
        # Also, let the user define a threshold for this! and take the maximum of the above and this value!

    # Get the additional information from the bunch object
    if re_deistope_precursor:
        
        precursor_info = bunch.precursor_information

        if not trust_precursor_scan_number:
            #Correct the precursor scan number to the correct one! This is sometimes incorrectly loaded from the raw file by ThermoRawLoader
            precursor_scan_number = bunch.source.find_previous_ms1(bunch.index)._data.scan_number
            precursor_scan_id = bunch.source.get_scan_by_index(precursor_scan_number).id        
            precursor_info.precursor_scan_id  = precursor_scan_id
            #print(precursor_info)

        _,re_deisotope_success, = bunch.precursor_information.find_monoisotopic_peak(
            trust_charge_state=False,
            #precursor_scan=bunch.source.find_previous_ms1(bunch.index),
            averageine=avgine,
            scorer=scorer,
            max_missed_peaks=max_missed_peaks,
            charge_range=charge_range,
            truncate_after=0.95,
            scale=scale,
            error_tolerance=error_tol,
            #error_tolerance=2e-6,
            minimum_intensity=0.,
            scale_method=scale_method,
        )

        #It can happen that the `find_monoisotopic_peak` is not complete/good enough
        #In that case, we also need a function that deisotopes the MS1 precurosor scan explicitly!
        #This is done by the following function: 
        # 
        # #TODO: Implement this function!
        def desitope_ms1_precusorsor():
            precursor_scan_id = bunch.precursor_information.id
            ms1_peakset = bunch.source.get_scan_by_id(precursor_scan_id).peak_set
            
            #Select the window of the precursor scan!
            ms1_range = ms1_peakset.between(bunch.isolation_window.lower_bound,bunch.isolation_window.higher_bound)
            
            ms1_ditp = ms_ditp.deconvolute_peaks(peaklist=ms1_range,
                                    averageine=avgine,
                                    scorer=scorer,
                                    max_missed_peaks=max_missed_peaks,
                                    charge_range=charge_range,
                                    truncate_after=truncate_after,
                                    scale=scale,
                                    error_tolerance=error_tol,
                                    minimum_intensity=minimum_intensity,
                                    scale_method=scale_method,
                                    )
            return ms1_ditp[0]
            #Extract the "additional information" from this MS1 scan further!

        add_info = re_deisotope_get_additional_bunch_info(bunch,bunch.precursor_information,re_deisotope_success)
    else:
        add_info = get_additional_bunch_info(bunch)

    deconvoluted_peaks = ms_ditp.deconvolute_peaks(
        bunch,
        averageine=avgine,
        scorer=scorer,
        max_missed_peaks=max_missed_peaks,
        charge_range=charge_range,
        truncate_after=truncate_after,
        scale=scale,
        error_tolerance=error_tol,
        minimum_intensity=minimum_intensity,
        scale_method=scale_method,
    )

    if return_data_frame:
        df = utils.convert_peaks_dataFrame(deconvoluted_peaks.peak_set,deisotoped=True)
        add_info_df = pl.concat([add_info] * len(deconvoluted_peaks.peak_set), how="vertical")
        df = pl.concat([df, add_info_df], how="horizontal")
        return df
    else:
        return deconvoluted_peaks, add_info


def deisotope_all_ms2_scans(
    raw_file_path,
    prune_single_charge_states=None,
    prune_single_charge_states_min_charge=1,
    prune_single_charge_states_keep_solo_unity_charge=True,
    prune_single_charge_states_keep_solo_unity_charge_mass_cutoff=900,
    aggregate_masses=None,
    aggregate_masses_index_column="neutral_mass",
    aggregate_masses_type="simple",
    aggregate_masses_tolerance=0.5,
    aggregate_masses_tolerance_type="absolute",
    deisotoping_parameters=dict(),
    re_deisotoping_parameters=dict(),
    min_num_peaks_per_ms2_scan=0,
):
    # Load the file
    raw_file_read = ms_ditp.data_source.thermo_raw_net.ThermoRawLoader(
        raw_file_path, _load_metadata=True
    )

    # List to store the deconvoluted results from each MS2 scan
    decon_list = []

    # Deconvolute each scan in the iteratable which is an MS2 scan!
    raw_file_read.make_iterator(grouped=False)
    bunch = next(raw_file_read)

    if re_deisotoping_parameters == {}:
        re_deisotoping_parameters = deisotoping_parameters

    while bunch is not None:
        if bunch.ms_level == 2:
            bunch.pick_peaks()

            # Deconvolute the next scan and get additional information about the precursors etc!
            t, add_info = deconvolute_ms2_bunch(
                bunch=bunch,
                parameters=deisotoping_parameters,
                re_deistope_precursor=True,
                re_deisotoping_parameters=re_deisotoping_parameters,
            )

            if len(t.peak_set) > min_num_peaks_per_ms2_scan:
                # Pruning masses which are present in only single charge states!
                if prune_single_charge_states is not None:
                    t2 = utils.keep_multiple_charge_states_only(
                        peakset=t.peak_set,
                        min_charge=prune_single_charge_states_min_charge,
                        keep_solo_unity_charge=prune_single_charge_states_keep_solo_unity_charge,
                        keep_solo_unity_charge_mass_cutoff=prune_single_charge_states_keep_solo_unity_charge_mass_cutoff,
                    )

                    # decon_list.append(utils.convert_peaks_dataFrame(t2))
                    add_info = pl.concat([add_info] * len(t2), how="vertical")
                    decon_list.append(
                        pl.concat(
                            [utils.convert_peaks_dataFrame(t2,deisotoped=True), add_info],
                            how="horizontal",
                        )
                    )
                    print(
                        "Scan number %d processed with %d deconvoluted peaks"
                        % (bunch.index, len(t2))
                    )

                else:
                    add_info = pl.concat([add_info] * len(t.peak_set), how="vertical")
                    decon_list.append(
                        pl.concat(
                            [utils.convert_peaks_dataFrame(t.peak_set,deisotoped=True), add_info],
                            how="horizontal",
                        )
                    )
                    print(
                        "Scan number %d processed with %d deconvoluted peaks"
                        % (bunch.index, len(t.peak_set))
                    )

        try:
            bunch = next(raw_file_read)
        except:
            break

    # Convert to a cumulative dataframe
    decon_df = pl.DataFrame()
    for i in decon_list:
        decon_df = decon_df.vstack(i)
    decon_df = decon_df.sort(pl.col("neutral_mass"))

    # Aggregate the masses based on some tolerance
    if aggregate_masses is not None:
        if aggregate_masses_type == "simple":
            decon_df = utils.custom_aggregate(
                decon_df,
                index_column=aggregate_masses_index_column,
                tolerance=aggregate_masses_tolerance,
                tolerance_type=aggregate_masses_tolerance_type,
            )
        elif aggregate_masses_type == "all":
            decon_df = utils.custom_aggregate_aggregate_all(
                decon_df,
                index_column=aggregate_masses_index_column,
                tolerance=aggregate_masses_tolerance,
                tolerance_type=aggregate_masses_tolerance_type,
            )

    return decon_df
