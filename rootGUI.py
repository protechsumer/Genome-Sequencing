from tkinter import *
from tkinter import messagebox
from tkinter import filedialog
import pandas as pd
from pandas import Series,DataFrame

from function import *

filename = " "
"""Functions"""
def openFile():
    global filename
    filename = filedialog.askopenfilename(initialdir="/", title="Select File",filetypes=(("text files",".txt"),("all files","*")))
    print(filename)

def mainProgram():

     start(filename)
     return

def fileOpen():
    print("Opening The File")

def aboutMe():
    about="Hello! I am Gene Analysis Interface."
    messagebox.showinfo("About Me",about)

def qut():
    quit()
""" Functions End"""

root = Tk()
root.title("Gene Analysis Interface")
root.geometry('700x500+250+100')
"""Menu"""
menu = Menu(root)
root.config(menu=menu)
subMenu = Menu(menu)
menu.add_cascade(label="Select File", menu=subMenu)
subMenu.add_command(label="Select New Fasta File", command=openFile)
subMenu.add_command(label="About Me", command=aboutMe)
subMenu.add_separator()
subMenu.add_command(label="Exit", command=qut)
"""MEnu End"""

"""Inside"""
topFrame = Frame(root)
label1 = Label(topFrame, text="Want to Select the Fasta File ->", bg="Red")
label2 = Label(topFrame, text="Want to start the process? ->", bg="Green")

button1 = Button(topFrame, text="Click Me", command=openFile)
button2 = Button(topFrame, text="Press Me", command=mainProgram)
label1.grid(row=1, column=1, sticky=E)
button1.grid(row=1, column=2, sticky=W)
label2.grid(row=2,column=1, sticky=E)
button2.grid(row=2, column=2, sticky=W)
topFrame.pack()
bottomFrame = Frame(root)
bottomFrame.pack(side=BOTTOM)
print(filename)
root.mainloop()