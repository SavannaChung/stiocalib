import numpy as np
import pandas as pd
import os
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows
import easygui as eg


# Read the excel file in dataframe
def read_data_meas_ctno(filename):
    """ read the meas ct no.
    -------------------------------input-------------------------------------
    filename (string) = excel file names
    insert_name : a string telling which excel sheet is read
    -------------------------------output-------------------------------------
    meas_ctno: {insert: meas HU, ...}
    """

    df = pd.concat(pd.read_excel(filename, sheet_name = None))
    df = df[["insert_phantom","inserts","hu_mean"]]

    # >>> get the categories
    phantoms = list(df["insert_phantom"].dropna().unique())

    print(f"phantoms: {phantoms, type(phantoms)}")

    values_lt = []
    index = []
    for phantom in phantoms:
        df1 = df[df.insert_phantom == phantom].groupby(["inserts"]).mean()

        iterables = [[phantom], list(df1.index)]
        index.append(pd.MultiIndex.from_product( iterables, names = ["phantoms", "inserts"]))
        values_lt.append(np.concatenate((df1.values), axis = None))

    index = np.concatenate(index)
    phantoms, inserts = zip(*index)
    values_lt = np.concatenate(values_lt)
    df1 = pd.DataFrame({"insert_phantom": phantoms, "insert": inserts, "hu_mean": values_lt})

    return df1


## >>> output the insert mean_hu from different phantom

def data_to_excel(df):

    name = "meas_ctno.xlsx"
    sheet_name = "meas_ctno"

    wb = Workbook()
    ws = wb.active
    ws.title = sheet_name

    for r in dataframe_to_rows(df, index = 0):
        ws.append(r)

    wb.save(name)


    # with pd.ExcelWriter(name) as excel_writer:
    #     df.to_excel(excel_writer, sheet_name = sheet_name)
    #     excel_writer.save()

    return


excel_file_path = eg.fileopenbox("Please select the excel file containing the measured HU data of inserts.")
path, file = os.path.split(excel_file_path)
os.chdir(path)

df = read_data_meas_ctno(file)
# print(df)
data_to_excel(df)

print(f"the measured ctno data have been written out.")
