import ms_deisotope as ms_ditp
import numpy as np
import polars as pl
import importlib


# Util functions

def convert_peaks_dataFrame(peak_set,deisotoped=False):
    """
    Converts a set of peak objects into a polar DataFrame.
    Parameters:
    peak_set (list): A list of ms_deisotope peak objects. Each object in the list should have attributes that contain data.
    Returns:
    pl.DataFrame: A polar DataFrame where each column corresponds to an attribute of the peak objects and each row corresponds to a peak object.
    """

    # Get the fields of the peak_set which contain data and are not callable or functions or init methods
    fields = [
        attr
        for attr in dir(peak_set[0])
        if not callable(getattr(peak_set[0], attr))
        and not attr.startswith("__")
        and not attr.startswith("_")
    ]

    data_frame = pl.DataFrame(
        {field: [getattr(peak, field) for peak in peak_set] for field in fields}
    )

    if deisotoped:
        # Cast the envelope to a string and the fit to its score and float vlaues
        envelope = [peak_set[i].envelope.pairs for i in range(len(peak_set))]
        envelope_str = [
            ";".join([f"({x[0]}, {x[1]})" for x in envelope[i]])
            for i in range(len(envelope))
        ]

        fit_score = [peak_set[i].fit.score for i in range(len(peak_set))]

        fit_npeaks = [peak_set[i].fit.npeaks for i in range(len(peak_set))]

        data_frame = data_frame.with_columns(
            pl.Series(name="envelope", values=envelope_str, strict=False),
            pl.Series(name="fit_score", values=fit_score, strict=False),
            pl.Series(name="fit_npeaks", values=fit_npeaks, strict=False),
        )
        if "index" in data_frame.columns:
            data_frame = data_frame.drop("index")
        if "fit" in data_frame.columns:
            data_frame = data_frame.drop("fit")
        if "chosen_for_msms" in data_frame.columns:
            data_frame = data_frame.drop("chosen_for_msms")

    return data_frame


def keep_multiple_charge_states_only(
    peakset,
    min_charge=1,
    keep_solo_unity_charge=True,
    keep_solo_unity_charge_mass_cutoff=np.inf,
    data_frame=False,
):
    """
    Filters a set of peaks to retain only those that appear at multiple charge states. Helps to increase confidence if the same peak is present at multiple charge states.
    Parameters:
    -----------
    peakset : iterable
        A collection of peak objects to be filtered.
    min_charge : int, optional
        Minimum number of charge states a peak must have to be retained. Default is 1.
    keep_solo_unity_charge : bool, optional
        If True, retains peaks with a single charge state of ±1 if their neutral mass is below a specified cutoff.
        Single nucleotiodes can only be present with a charge state of 1.
        So below a certai neutral mass value, we accept it even if only a single charge state with charge == +-1 is present.
        Default is True.
    keep_solo_unity_charge_mass_cutoff : float, optional
        The neutral mass cutoff below which peaks with a single charge state of ±1 are retained. Default is np.inf.
    data_frame : bool, optional
        If True, returns the filtered peaks as a DataFrame. Otherwise, returns a set of peaks. Default is False.
    Returns:
    --------
    set or DataFrame
        A set of peaks that meet the filtering criteria, or a DataFrame if `data_frame` is True.

    #TODO: Change cutoff to some factor*heaviest nucleoside mass!
    """

    seen_peaks = set()
    multi_charge_peaks = set()

    for peak in peakset:
        if peak in seen_peaks:
            continue

        seen_peaks.add(peak)

        all_for = peakset.all_peaks_for(peak.neutral_mass)

        if (
            len(all_for) > min_charge
        ):  # Only store if the same peak is present at multiple m/z values, i.e. a proxy for multiple charge states
            all_for = sorted(all_for, key=lambda x: x.mz, reverse=1)

            for p in all_for:
                seen_peaks.add(p)
                multi_charge_peaks.add(p)

        elif keep_solo_unity_charge:
            if len(all_for) == 1:
                if (
                    abs(all_for[0].charge) == 1
                    and all_for[0].neutral_mass < keep_solo_unity_charge_mass_cutoff
                ):
                    for p in all_for:
                        seen_peaks.add(p)
                        multi_charge_peaks.add(p)

    if data_frame:
        return convert_peaks_dataFrame(list(multi_charge_peaks))
    else:
        return multi_charge_peaks


# Custom aggregation function
def custom_aggregate(
    df,
    tolerance=0.1,
    index_column="neutral_mass",
    name_intesity_col="intensity",
    tolerance_type="absolute",
    additional_matching_criterion=None,
):
    """
    Aggregates rows in a DataFrame based on a specified tolerance and optional additional matching criterion.
    Parameters:
    df (DataFrame): The input DataFrame to aggregate.
    tolerance (float, optional): The tolerance value for grouping rows. Default is 0.1.
    index_column (str, optional): The column name to use for grouping. Default is "neutral_mass".
    name_intesity_col (str, optional): The column name for intensity values. Default is "intensity".
    tolerance_type (str, optional): The type of tolerance to use ('absolute' or 'relative'). Default is 'absolute'.
    additional_matching_criterion (str, optional): An additional column name to use for matching rows. Default is None.
    Returns:
    DataFrame: A new DataFrame with aggregated rows based on the specified tolerance and additional matching criterion.

    """

    def deviation_check(value, ref_value):
        if tolerance_type == "absolute":
            if abs(value - ref_value) <= tolerance:
                return True
            else:
                return False
        elif tolerance_type == "relative":
            if abs(value / ref_value - 1) <= tolerance:
                return True
            else:
                return False

    df = df.sort(index_column)
    index_column_number = df.columns.index(index_column)

    if additional_matching_criterion is not None:
        additional_column_number = df.columns.index(additional_matching_criterion)

    groups = []
    current_group = []
    for i, row in enumerate(df.rows()):
        if not current_group:
            current_group.append(row)
        else:
            if additional_matching_criterion is not None:
                if (
                    deviation_check(
                        row[index_column_number], current_group[0][index_column_number]
                    )
                    and row[additional_column_number]
                    == current_group[0][additional_column_number]
                ):
                    current_group.append(row)
                else:
                    groups.append(current_group)
                    current_group = [row]
            else:
                if deviation_check(
                    row[index_column_number], current_group[0][index_column_number]
                ):
                    current_group.append(row)
                else:
                    groups.append(current_group)
                    current_group = [row]

    if current_group:
        groups.append(current_group)

    aggregated_data = []

    for group in groups:
        # index_column_avg = sum(row[index_column_number]*] for row in group) / len(group)
        intensity_sum = sum(row[df.columns.index(name_intesity_col)] for row in group)
        index_column_avg = (
            sum(
                row[index_column_number] * row[df.columns.index(name_intesity_col)]
                for row in group
            )
            / intensity_sum
        )

        if additional_matching_criterion is not None:
            aggregated_data.append(
                (index_column_avg, intensity_sum, group[0][additional_column_number])
            )
        else:
            aggregated_data.append((index_column_avg, intensity_sum))

    if additional_matching_criterion is not None:
        return pl.DataFrame(
            aggregated_data,
            schema=[index_column, name_intesity_col, additional_matching_criterion],
            orient="row",
        )
    else:
        return pl.DataFrame(
            aggregated_data, schema=[index_column, name_intesity_col], orient="row"
        )


def custom_aggregate_aggregate_all(
    df,
    tolerance,
    index_column="neutral_mass",
    name_intesity_col="intensity",
    tolerance_type="absolute",
    additional_matching_criterion=None,
):
    """
    Aggregates rows in a DataFrame based on a specified tolerance and optional additional matching criterion. Uses polars aggregation to aggregate other properties and columns as well.
    Parameters:
    df (DataFrame): The input DataFrame to aggregate.
    tolerance (float, optional): The tolerance value for grouping rows. Default is 0.1.
    index_column (str, optional): The column name to use for grouping. Default is "neutral_mass".
    name_intesity_col (str, optional): The column name for intensity values. Default is "intensity".
    tolerance_type (str, optional): The type of tolerance to use ('absolute' or 'relative'). Default is 'absolute'.
    additional_matching_criterion (str, optional): An additional column name to use for matching rows. Default is None.
    Returns:
    DataFrame: A new DataFrame with aggregated rows based on the specified tolerance and additional matching criterion.
    """

    def deviation_check(value, ref_value):
        if tolerance_type == "absolute":
            if abs(value - ref_value) <= tolerance:
                return True
            else:
                return False
        elif tolerance_type == "relative":
            if abs(value / ref_value - 1) <= tolerance:
                return True
            else:
                return False

    df = df.sort(index_column)
    index_column_number = df.columns.index(index_column)

    if additional_matching_criterion is not None:
        additional_column_number = df.columns.index(additional_matching_criterion)

    groups = []
    current_group = []
    index_list = [0] * len(df)
    index = 0

    for i, row in enumerate(df.rows()):
        if not current_group:
            current_group.append(row)
            index_list[i] = index
        else:
            if additional_matching_criterion is not None:
                if (
                    deviation_check(
                        row[index_column_number], current_group[0][index_column_number]
                    )
                    and row[additional_column_number]
                    == current_group[0][additional_column_number]
                ):
                    current_group.append(row)
                    index_list[i] = index
                else:
                    groups.append(current_group)
                    current_group = [row]
                    index += 1
                    index_list[i] = index
            else:
                if deviation_check(
                    row[index_column_number], current_group[0][index_column_number]
                ):
                    current_group.append(row)
                    index_list[i] = index
                else:
                    groups.append(current_group)
                    current_group = [row]
                    index += 1
                    index_list[i] = index

    if current_group:
        groups.append(current_group)

    df = df.with_columns(pl.Series(name="group_index", values=index_list))

    col_names = [
        str(col)
        for col in df.columns
        if col not in [index_column, name_intesity_col, "group_index"]
    ]
    # col_names_dtype_objects = [col for col in col_names if df.schema[col] == type(df.schema["fit"])]

    agg_df = (
        df.group_by("group_index")
        .agg(  # [pl.mean(index_column)] +
            [
                (pl.col(index_column) * pl.col(name_intesity_col)).sum()
                / pl.col(name_intesity_col).sum()
            ]
            + [pl.col(name_intesity_col).alias("individual_intensities")]
            + [pl.sum(name_intesity_col)]
            + [pl.col(col) for col in col_names]
            # [pl.col(col) for col in col_names if col not in col_names_dtype_objects]
        )
        .sort(pl.col(index_column))
    )

    return agg_df
