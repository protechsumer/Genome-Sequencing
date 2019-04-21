from tkinter import messagebox
from tkinter import filedialog
import pandas as pd
import numpy as np
import array as arr
import csv
from pandas import Series,DataFrame
from collections import Counter

import mysql.connector
mydb=mysql.connector.connect(
    host="localhost",
    user="root",
    password="root",
    database="gs"
)
mycursor=mydb.cursor()


def errorchecking(file):
    k = 0
    l = 0
    for line in file:
        if (line[0] == '>'):
            if (line[len(line) - 1] == '\n'):
                k = k + 1
            else:
                k = -1
                break
        else:
            if (line[0] == '\n'):
                k = -1
                break
            for i in line:
                if (i != '\n'):
                    if ((i != 'A') & (i != 'T') & (i != 'G') & (i != 'C')):
                        l = -1
                        break
    #file.close()
    #if ((k == -1) | (l == -1)):
    #    return 0
    #return 1
    if ((k == -1) | (l == -1)):
        errmsg = "File Has Error. Select another file"
        print(errmsg)
        messagebox.showinfo("Error", errmsg)
        file.close()
        return 0
    else:
        accpmsg = "File is Correct."
        messagebox.showinfo("Correct", accpmsg)
        file.close()
        return 1


def dataframeform(file):
    k = errorchecking(file)
    if (k == 0):
        return 0
    print("File Ok")
    file = open('ecoli.txt', 'r')
    df = pd.DataFrame(columns=['Slno', 'Information', 'Gene'])
    arr = ""
    info = ""
    k = -1
    for line in file:
        if (line[0] == '>'):
            k += 1
            df = df.append({'Slno': k, 'Information': info, 'Gene': arr}, ignore_index=True)
            info = line[1:len(line) - 1]
            arr = ""
        else:
            arr = arr + line[0:len(line) - 1]
    df = df.append({'Slno': k + 1, 'Information': info, 'Gene': arr}, ignore_index=True)
    df = df.loc[1:]
    return df


def atgccount(df):
    grr = df.Gene
    da = pd.DataFrame(columns=['A', 'C', 'G', 'T', 'Length', '(G+C)%'])
    da = da.append({'A': 0, 'C': 0, 'G': 0, 'T': 0, 'Length': 0, '(G+C)%': 0}, ignore_index=True)
    for i in range(1, len(grr) + 1):
        a = list(grr[i])
        uniq, counts = np.unique(a, return_counts=True)
        da = da.append({'A': counts[0], 'C': counts[1], 'G': counts[2], 'T': counts[3], 'Length': sum(counts),
                        '(G+C)%': ((counts[1] + counts[2]) / sum(counts)) * 100}, ignore_index=True)
    db = pd.concat([df, da], axis=1, sort=False)
    db = db.loc[1:]
    return db


def eliminate(arr, a, b):
    l = 0
    for i in range(a, b):
        if (arr[i] == ','):
            l = 1
            break
    return l


def elimination(db):
    df = pd.DataFrame(db)
    crr = list(db['Information'])
    p = 0
    for i in range(0, len(crr)):
        arr = crr[i]
        f = 0
        for j in range(0, len(arr)):
            if (arr[j] == ':'):
                for k in range((j + 1), len(arr)):
                    if (arr[k] == ' '):
                        f = eliminate(arr, j + 1, k)
                        if (f == 1):
                            df = df.drop(df.index[i - p])
                            p = p + 1
                        k = len(arr) + 1
                        j = len(arr) + 1
    return df


def binarysearch(arr, l, r, x):
    while (l <= r):
        m = l + (r - l) // 2
        if (arr[m] == x):
            return m
        if (arr[m] < x):
            l = m + 1
        else:
            r = m - 1
    return -1


def rangefinder(df):
    df = elimination(df)
    data = pd.read_csv("genelist.csv")
    db = pd.DataFrame(columns=['Slno', 'Information', 'Gene_Sequence', 'A', 'C', 'G', 'T', 'Length(ATGC)', '(G+C)%',
                               'Location1', 'Location2', 'Strand', 'Length', 'PID', 'Gene', 'Synonym',
                               'Code', 'COG', 'Product'])
    info = list(df['Information'])
    gene = list(df['Gene'])
    A = list(df['A'])
    C = list(df['C'])
    G = list(df['G'])
    T = list(df['T'])
    Length = list(df['Length'])
    GC = list(df['(G+C)%'])
    p = 0
    for i in range(0, len(info)):
        arr = info[i]
        prr = ""
        qrr = ""
        b = 0
        for j in range(0, len(arr)):
            if (arr[j] == ':'):
                n = 0
                m = 0
                c = 0
                if (arr[j + 1] == 'c'):
                    c = 1
                if (c == 0):
                    for k in range((j + 1), len(arr)):
                        if (arr[k] == ' '):
                            break
                        if (arr[k] == '-'):
                            n = 1
                            m = 0
                        if ((n == 0) & (arr[k] != '-')):
                            prr = prr + arr[k]
                            m = m + 1
                        if ((n == 1) & (arr[k] != '-')):
                            qrr = qrr + arr[k]
                            m = m + 1
                if (c == 1):
                    for k in range((j + 2), len(arr)):
                        if (arr[k] == ' '):
                            break
                        if (arr[k] == '-'):
                            n = 1
                            m = 0
                        if ((n == 0) & (arr[k] != '-')):
                            qrr = qrr + arr[k]
                            m = m + 1
                        if ((n == 1) & (arr[k] != '-')):
                            prr = prr + arr[k]
                            m = m + 1
        prr = (int(prr))
        qrr = (int(qrr))
        b = binarysearch(data['Location1'], 0, len(data['Location1']), prr)
        if (b != -1):
            if (qrr == data['Location2'][b]):
                p = p + 1
                db = db.append(
                    {'Slno': p, 'Information': info[i], 'Gene_Sequence': gene[i], 'A': A[i], 'C': C[i], 'G': G[i],
                     'T': T[i], 'Length(ATGC)': Length[i], '(G+C)%': GC[i], 'Location1': data['Location1'][b],
                     'Location2': data['Location2'][b], 'Strand': data['Strand'][b], 'Length': data['Length'][b],
                     'PID': data['PID'][b], 'Gene': data['Gene'][b], 'Synonym': data['Synonym'][b],
                     'Code': data['Code'][b], 'COG': data['COG'][b], 'Product': data['Product'][b]}, ignore_index=True)
            else:
                if (((b - 1) >= 0) & ((b + 1) < len(info))):
                    if (qrr == data['Location2'][b - 1]):
                        db = db.append(
                            {'Slno': p, 'Information': info[i], 'Gene_Sequence': gene[i], 'A': A[i], 'C': C[i],
                             'G': G[i], 'T': T[i], 'Length(ATGC)': Length[i], '(G+C)%': GC[i],
                             'Location1': data['Location1'][b - 1], 'Location2': data['Location2'][b - 1],
                             'Strand': data['Strand'][b - 1], 'Length': data['Length'][b - 1],
                             'PID': data['PID'][b - 1], 'Gene': data['Gene'][b - 1], 'Synonym': data['Synonym'][b - 1],
                             'Code': data['Code'][b - 1], 'COG': data['COG'][b - 1], 'Product': data['Product'][b - 1]},
                            ignore_index=True)
                        p = p + 1
                    if (qrr == data['Location2'][b + 1]):
                        db = db.append(
                            {'Slno': p, 'Information': info[i], 'Gene_Sequence': gene[i], 'A': A[i], 'C': C[i],
                             'G': G[i], 'T': T[i], 'Length(ATGC)': Length[i], '(G+C)%': GC[i],
                             'Location1': data['Location1'][b + 1], 'Location2': data['Location2'][b + 1],
                             'Strand': data['Strand'][b + 1], 'Length': data['Length'][b + 1],
                             'PID': data['PID'][b + 1], 'Gene': data['Gene'][b + 1], 'Synonym': data['Synonym'][b + 1],
                             'Code': data['Code'][b + 1], 'COG': data['COG'][b + 1], 'Product': data['Product'][b + 1]},
                            ignore_index=True)
                        p = p + 1
    return db


def writingfile(db):
    #print("Do you want to make the information present in dataframe and .csv file ? : ")
    #n=input()
    filename2=filedialog.asksaveasfile(initialdir="/",title="Save File As",filetypes=(("csv files","*.csv"),("all files",".*")))
    print(filename2)
    #if(n.lower()=='y'):
    #  filename=input()
    #  filename=filename+".csv"

    db.to_csv(filename2)

    for index, row in db.iterrows():
        mycursor.execute(
            "INSERT INTO genome(SL_No,Information,Gene_Sequence,A_Count,T_Count,G_Count,C_Count,Gene_Length,GC_Percent,Location1,Location2,Strand,Length,PID,Gene,Synonym,Code,COG,Product) \
                values(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s) ON DUPLICATE KEY UPDATE SL_No=%s", \
            (row['Slno'], row['Information'], row['Gene_Sequence'], row['A'], row['T'], row['G'], row['C'], row['Length(ATGC)'], row['(G+C)%'],row['Location1'],row['Location2'],\
            row['Strand'],row['Length'],row['PID'],row['Gene'],row['Synonym'],row['Code'],row['COG'],row['Product'],row['Slno']))
        mydb.commit()

    # print("The file has been written with name",filename)

    mycursor.execute("SELECT * FROM genome")
    print("PRINTING FROM DATABASE::")
    for i in mycursor:
        print(i)



def start(filename):
    file = open(filename, 'r')
    df = dataframeform(file)
    if (type(df) == int):
        if (df == 0):
            print("File error")
            return 0
    # print(df)
    df = atgccount(df)
    # print(df)
    df = rangefinder(df)
    # print(df)
    writingfile(df)


#start(filename)