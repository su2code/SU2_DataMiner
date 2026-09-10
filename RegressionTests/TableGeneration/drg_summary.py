"""Summarize a written Dragon table so its values, and not only its structure, are regression tested.

The regression harness splits each line on commas and compares the resulting words as floats. Data
lines of a Dragon table are tab separated, so they reach the comparison as a single non-numeric word
and are never compared. Reducing the table to one comma separated line per variable puts the values
themselves under the same tolerance as the rest of the suite.
"""
import numpy as np


def readTableData(drg_path:str):
    """Read the variable names and the table data of a Dragon table file.

    :param drg_path: path of the Dragon table file.
    :type drg_path: str
    :return: variable names and the table data of every level stacked together.
    :rtype: tuple[list[str], np.ndarray[float]]
    """
    with open(drg_path, "r") as fid:
        lines = [line.strip() for line in fid]

    n_variables = int(lines[lines.index("[Number of variables]") + 1])
    names_start = lines.index("[Variable names]") + 1
    variable_names = [lines[i].split(":", 1)[-1] for i in range(names_start, names_start + n_variables)]

    table_data = []
    for line in lines[lines.index("<Data>") + 1:lines.index("</Data>")]:
        if line and not line.startswith("<"):
            table_data.append([float(value) for value in line.split()])

    return variable_names, np.asarray(table_data, dtype=float)


def writeTableSummary(drg_path:str, summary_path:str):
    """Write the extent of every variable of a Dragon table as a comma separated summary.

    :param drg_path: path of the Dragon table file to summarize.
    :type drg_path: str
    :param summary_path: path of the summary file to write.
    :type summary_path: str
    """
    variable_names, table_data = readTableData(drg_path)

    with open(summary_path, "w") as fid:
        fid.write("variable,minimum,maximum,mean\n")
        for iVar, name in enumerate(variable_names):
            fid.write("%s,%+.10e,%+.10e,%+.10e\n" % (name, np.min(table_data[:, iVar]), \
                                                     np.max(table_data[:, iVar]), \
                                                     np.mean(table_data[:, iVar])))
    return
