import numpy as np
import scipy
import scipy.linalg as la
from tkinter import *
import tkinter.font as font
import pprint
import re
from functions import *

class Gui:
    method=None
    numberOfVar=None
    b=[]
    xu=None
    xl=None
    x0=0
    x1=0
    initial_Guessstr=None
    initial_Guess=[]
    equationStr=None
    matrix= []
    function=None
    stopping_condition=None
    noOfitrations=50
    tolerance= 0.00001
    def get_initialGuess(self):
        self.initial_GuessStr=self.initial_GuessStr.split()
        for x in self.initial_GuessStr:
            self.initial_Guess.append(float(x))
        return

    def get_coeff(self):
        arrayeqs=self.equationStr.split("\n")
        for i in range(self.numberOfVar):
            self.b.append(float(arrayeqs[i][arrayeqs[i].find("=")+1:]))
            mystring = arrayeqs[i][:arrayeqs[i].find("=")]
            regex_coeff = '([+-]\d*\.{0,1}\d+) x'
            coeffs = [float(arrayeqs[i][:arrayeqs[i].find(" ")])] + [float(x) for x in re.findall(regex_coeff, mystring)]
            self.matrix.append(coeffs)
        print("ffffffffffffff")
        print(self.matrix,self.b)
        return


    def solve_btn(self):
        gui_sol=Tk()
        gui_sol.geometry("400x400")
        if self.method=="Gauss Elimination":
            x, timex =GEPP(self.matrix,self.b,self.numberOfVar)
            Label(gui_sol, text="x=" + str(x)).place(x=25, y=25)
            Label(gui_sol, text="time: " + str(timex)).place(x=25, y=50)

        elif self.method=="Gauss-jordan":
            aug=[]
            for i in range(self.numberOfVar):
                aug.append(self.matrix[i])
                aug[i].append(self.b[i])
            x, timex = gaussJordan(aug, self.matrix, self.numberOfVar)
            Label(gui_sol, text="x=" + str(x)).place(x=25, y=25)
            Label(gui_sol, text="time: " + str(timex)).place(x=25, y=50)
        elif self.method=="LU Decomposition":
            x, timex = luDecomposition(self.matrix,self.b,self.numberOfVar)
            if x==0:
                Label(gui_sol, text="no solution").place(x=25, y=25)
            elif x==2:
                Label(gui_sol, text="Infimit number of solutions").place(x=25, y=25)
            else:
                Label(gui_sol, text="x=" + str(x)).place(x=25, y=25)
                Label(gui_sol, text="time: " + str(timex)).place(x=25, y=50)
        elif self.method=="Gauss-seidil":
            x,nu,timex,error=gauss_Seidil(self.matrix,self.b,self.initial_Guess,self.tolerance,self.noOfitrations)
            Label(gui_sol, text="x=" + str(x)).place(x=25, y=25)
            Label(gui_sol, text="time: " + str(timex)).place(x=25, y=50)
            Label(gui_sol, text="# itration: " + str(nu)).place(x=25, y=75)
            Label(gui_sol, text="error: " + str(error)).place(x=25, y=100)
        elif self.method=="jacobi-iteration":
            print(self.matrix,self.b, self.initial_Guess, self.tolerance, self.noOfitrations)
            x,nu,timex,error=jacobi_iterations(self.matrix,self.b, self.initial_Guess, self.tolerance, self.noOfitrations)
            Label(gui_sol, text="x=" + str(x)).place(x=25, y=25)
            Label(gui_sol, text="time: " + str(timex)).place(x=25, y=50)
            Label(gui_sol, text="# itration: " + str(nu)).place(x=25, y=75)
            Label(gui_sol, text="error: " + str(error)).place(x=25, y=100)
        elif self.method=="bisection":
            xr,f_xr,nu,timex,error=bracketing(1,self.function,self.xl,self.xu,self.tolerance,self.noOfitrations)
            print(bracketing(1,self.function,self.xl,self.xu,self.tolerance,self.noOfitrations))
            Label(gui_sol, text="x=" + str(xr)).place(x=25, y=25)
            Label(gui_sol, text="f(xr): " + str(timex)).place(x=25, y=50)
            Label(gui_sol, text="time: " + str(timex)).place(x=25, y=50)
            Label(gui_sol, text="# itration: " + str(nu)).place(x=25, y=75)
            Label(gui_sol, text="error: " + str(error)).place(x=25, y=100)
            plot_function(self.function, True)
        elif self.method=="False-position":
            xr,f_xr,nu,timex,error=bracketing(2, self.function, self.xl, self.xu, self.tolerance, self.noOfitrations)
            Label(gui_sol, text="x=" + str(xr)).place(x=25, y=25)
            Label(gui_sol, text="f(xr): " + str(timex)).place(x=25, y=50)
            Label(gui_sol, text="time: " + str(timex)).place(x=25, y=50)
            Label(gui_sol, text="# itration: " + str(nu)).place(x=25, y=75)
            Label(gui_sol, text="error: " + str(error)).place(x=25, y=100)
            plot_function(self.function, True)
        elif self.method=="Fixed-point":
            xr,nu,timex,error=fixedPointIteration(self.function,self.x0,self.tolerance,self.noOfitrations)
            Label(gui_sol, text="x=" + str(xr)).place(x=25, y=25)
            Label(gui_sol, text="f(xr): " + str(timex)).place(x=25, y=50)
            Label(gui_sol, text="time: " + str(timex)).place(x=25, y=50)
            Label(gui_sol, text="# itration: " + str(nu)).place(x=25, y=75)
            Label(gui_sol, text="error: " + str(error)).place(x=25, y=100)
            plot_function(self.function, False)
        elif self.method=="Newton-Raphson":
            x,timex=newtonRaphson(self.function,self.x0,self.tolerance,self.noOfitrations)
            Label(gui_sol, text="x=" + str(x)).place(x=25, y=25)
            Label(gui_sol, text="time: " + str(timex)).place(x=25, y=50)
            plot_function(self.function, True)
        elif self.method=="secant Method":
            x, timex = secant(self.function, self.x0,self.x1 ,self.tolerance, self.noOfitrations)
            Label(gui_sol, text="x=" + str(x)).place(x=25, y=25)
            Label(gui_sol, text="time: " + str(timex)).place(x=25, y=50)
            plot_function(self.function, True)

        gui_sol.mainloop()




root=Tk()
root.title("calculator")
root.geometry("400x400")
root.config(bg='#49A')

solve=Gui

methodList = [
    "Gauss Elimination",
    "Gauss-jordan",
    "LU Decomposition",
    "Gauss-seidil",
    "jacobi-iteration",
    "bisection",
    "False-position",
    "Fixed-point",
    "Newton-Raphson",
    "secant Method"
]
def select(x):
    method = clicked.get()
    solve.method = clicked.get()
    def functionEntry():
        global functionE
        if(solve.method=="Fixed-point"):
            txt1 = StringVar(root, value="Enter g(x)")
        else:
            txt1 = StringVar(root, value="Enter the function(with x var)")

        functionE=Entry(root,textvariable=txt1,width=35)
        functionE.place(x=75,y=50)

    def iteration_or_tolerance():
        global noOfiItEntry
        noOfiItEntry = Entry(root)
        noOfiItEntry.place(x=25, y=85)
        def stopping_cond(x):
            solve.stopping_condition = click.get()
            if (solve.stopping_condition == "Number of iterations"):
                noOfiItEntry.insert(0, "Enter # of itration")
            elif (solve.stopping_condition == "Abslute relative Error"):
                noOfiItEntry.insert(0, "enter abs relative error")
        click = StringVar(root)
        click.set("stop condition")
        Stop_condition = OptionMenu(root, click, "Number of iterations", "Abslute relative Error",command=stopping_cond)
        Stop_condition.config(width=20)
        Stop_condition.config(bg='light blue')
        Stop_condition.place(x=200, y=85)



    def get_paramerters():
        if(solve.method=="Gauss Elimination" or solve.method=="Gauss-jordan" or solve.method=="LU Decomposition" or solve.method=="Gauss-seidil"or solve.method=="jacobi-iteration"):
            nOfvar = int(noOfvarE.get())
            solve.numberOfVar=nOfvar
            solve.equationStr=equations.get("1.0", "end-1c")
            solve.get_coeff(solve)
        if(solve.method=="Gauss-seidil"or solve.method=="jacobi-iteration"):
            solve.initial_GuessStr=initialGuess.get()
            if(solve.stopping_condition=="Number of iterations"):
                solve.noOfitrations=int(noOfiItEntry.get())
            if (solve.stopping_condition == "Abslute relative Error"):
                solve.tolerance = float(noOfiItEntry.get())
            solve.get_initialGuess(solve)

        elif(solve.method=="bisection" or solve.method=="False-position"):
            solve.function=functionE.get()
            solve.xl=float(xlE.get())
            solve.xu=float(xuE.get())
            if (solve.stopping_condition == "Number of iterations"):
                solve.noOfitrations = noOfiItEntry.get()
            if (solve.stopping_condition == "Abslute relative Error"):
                solve.tolerance = noOfiItEntry.get()
        elif(solve.method=="Fixed-point" or solve.method=="Newton-Raphson" or solve.method=="secant Method"):
            solve.function = functionE.get()
            solve.x0=float(x0E.get())
            if (solve.stopping_condition == "Number of iterations"):
                solve.noOfitrations = float(noOfiItEntry.get())
            if (solve.stopping_condition == "Abslute relative Error"):
                solve.tolerance = float(noOfiItEntry.get())
            solve.x0=float(x0E.get())
        if(solve.method == "secant Method"):
            solve.x1=float(x1E.get())

        solve.solve_btn(solve)






    def getXlXU():
        global xlE
        global xuE
        textxl,textxu=StringVar(root,value="Enter Xl"),StringVar(root,value="Enter XU")
        xlE=Entry(root,textvariable=textxl)
        xlE.place(x=25,y=125)
        xuE=Entry(root,textvariable=textxu)
        xuE.place(x=200,y=125)

    def get_x0():
        global x0E
        textx0=StringVar(root,value="Enter x0")
        x0E=Entry(root,textvariable=textx0)
        x0E.place(x=130,y=130)
    def get_x1():
        global x1E
        textx1=StringVar(root,value="Enter x1")
        x1E=Entry(root,textvariable=textx1)
        x1E.place(x=130,y=150)



    if (method== "Gauss Elimination" or method=="Gauss-jordan" or method=="LU Decomposition" or method=="Gauss-seidil"or method=="jacobi-iteration"):
        txt= StringVar(root, value="enter number of variables")
        noOfvarE=Entry(root,textvariable=txt,width=25)
        noOfvarE.place(x=50,y=50)
        equations = Text(root, width=35, height=7)
        equations.place(x=50, y=120)


    if(method=="Gauss-seidil"or method=="jacobi-iteration"):
        txt5=StringVar(root, value="initial guess ")
        initialGuess=Entry(root,textvariable=txt5)
        initialGuess.place(x=225,y=50)
        iteration_or_tolerance()


    elif(method=="bisection" or method=="False-position"):
        functionEntry()
        iteration_or_tolerance()
        getXlXU()


    elif (method == "Fixed-point" or method == "Newton-Raphson" ):
        functionEntry()
        iteration_or_tolerance()
        get_x0()
    elif(method=="secant Method"):
        functionEntry()
        iteration_or_tolerance()
        get_x0()
        get_x1()


    f = font.Font(size=12)
    btn_solve = Button(root, text="solve", padx=50, pady=30, command=lambda: get_paramerters())
    btn_solve.place(x=100, y=250)
    btn_solve["font"] = f



clicked=StringVar()
clicked.set("choose method")
dropdown =OptionMenu(root,clicked,*methodList,command=select)
dropdown.config(width=40)
dropdown.config(bg='#191970')
dropdown.place(x=50,y=0)



root.mainloop()
