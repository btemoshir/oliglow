import ms_deisotope as ms_ditp


# Defines the averagine for RNA with and without backbone
def create_averageine(with_backbone=True,with_thiophosphate_backbone=False):
    """
    Calculate the average elemental composition of RNA nucleosides, optionally including the backbone.
    Args:
        with_backbone (bool): If True, include the phosphate group in the calculation. Defaults to True.
    Returns:
        dict: A dictionary representing the average elemental composition of RNA nucleosides.
              The keys are the elements ("C", "H", "N", "O", "P") and the values are the average counts of each element.
              If with_backbone is True, the composition includes the phosphate group.
    """

    if with_thiophosphate_backbone:
        # Nucleosides:
        Adenosine = {"C": 10, "H": 13, "N": 5, "O": 4, "P": 0, "S": 0}
        Cytidine = {"C": 9, "H": 13, "N": 3, "O": 5, "P": 0, "S": 0}
        Guanosine = {"C": 10, "H": 13, "N": 5, "O": 5, "P": 0, "S": 0}
        Uridine = {"C": 9, "H": 12, "N": 2, "O": 6, "P": 0, "S": 0}

        elements = ["C", "H", "N", "O", "P", "S"]
    
    else:
        # Nucleosides:
        Adenosine = {"C": 10, "H": 13, "N": 5, "O": 4, "P": 0}
        Cytidine = {"C": 9, "H": 13, "N": 3, "O": 5, "P": 0}
        Guanosine = {"C": 10, "H": 13, "N": 5, "O": 5, "P": 0}
        Uridine = {"C": 9, "H": 12, "N": 2, "O": 6, "P": 0}

        elements = ["C", "H", "N", "O", "P"]

    nucleosides_composition = [Adenosine, Cytidine, Guanosine, Uridine]

    average_nucleoside_rna = {key: 0 for key in elements}

    for e in elements:
        average_nucleoside_rna[e] = sum(
            [float(nucleosides[e]) for nucleosides in nucleosides_composition]
        ) / len(nucleosides_composition)
    # print("RNA nucleoside averagine = ", average_nucleoside_rna)

    # for e in elements:
    averagine_rna_with_backbone = average_nucleoside_rna.copy()
    averagine_rna_with_backbone["O"] += 2  # Add 2 oxygens for the phosphate group
    averagine_rna_with_backbone["P"] += 1  # Add 1 phosphorus for the phosphate group
    
    if with_thiophosphate_backbone:
        averagine_rna_with_thiophosphate_backbone = average_nucleoside_rna.copy()
        averagine_rna_with_thiophosphate_backbone["O"] += 1  # Add 1 oxygen for the phosphate group
        averagine_rna_with_thiophosphate_backbone["S"] += 1  # Add 1 sulphur for the thiophosphate group
        averagine_rna_with_thiophosphate_backbone["P"] += 1  # Add 1 phosphorus for the phosphate group

    if with_backbone:
        return averagine_rna_with_backbone
    elif with_thiophosphate_backbone:
        return averagine_rna_with_thiophosphate_backbone
    else:
        return average_nucleoside_rna


average_nucleoside_rna = ms_ditp.Averagine(create_averageine(with_backbone=False))
averagine_rna_with_backbone = ms_ditp.Averagine(create_averageine(with_backbone=True))
averagine_rna_with_thiophosphate_backbone = ms_ditp.Averagine(create_averageine(with_backbone=False,with_thiophosphate_backbone=True))
