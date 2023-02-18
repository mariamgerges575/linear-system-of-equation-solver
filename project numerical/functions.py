import time
from numpy import array , zeros, dot, diag
import re
import matplotlib.pyplot as plt
import sympy as smp
import datetime
import matplotlib.pyplot as plt
import numpy as np


def luDecomposition(mat, B, n):
    start = time.time()

    aug = np.zeros((n, n + 1))

    for i in range(n):
        for j in range(n + 1):
            if (j != n):
                aug[i][j] = mat[i][j]
            else:
                aug[i][j] = B[i]
    print("aug = \n", aug)

    # RANK OF MATRCIS
    rankA = np.linalg.matrix_rank(mat)
    rankAug = np.linalg.matrix_rank(aug)
    print("rankA = ", rankA, "\nrankAug = ", rankAug, "\n")
    if (rankA == rankAug):
        if (rankA == n):
            print("The system has a Unique solution")
            flag = 1
        else:
            print("The system has a infinte number of solution")
            flag = 2
    else:
        print("The system doesnot has any solution")
        flag = 0

    lower = [[0 for x in range(n)]
             for y in range(n)]
    upper = [[0 for x in range(n)]
             for y in range(n)]

    # Decomposing matrix into Upper
    # and Lower triangular matrix
    for i in range(n):

        # Upper Triangular
        for k in range(i, n):

            # Summation of L(i, j) * U(j, k)
            sum = 0
            for j in range(i):
                sum += (lower[i][j] * upper[j][k])

            # Evaluating U(i, k)
            upper[i][k] = mat[i][k] - sum
        # print(upper)

        # Lower Triangular
        for k in range(i, n):
            if (i == k):
                lower[i][i] = 1  # Diagonal as 1
            else:

                # Summation of L(k, j) * U(j, i)
                sum = 0
                for j in range(i):
                    sum += (lower[k][j] * upper[j][i])

                # Evaluating L(k, i)
                lower[k][i] = float((mat[k][i] - sum) / upper[i][i])

    print("lower triangular matrix :\n", lower)
    print("upper triangular matrix :\n", upper)

    lower = array(lower)
    upper = array(upper)

    print("l y = B\ny = \n")
    y = GEPPLU(lower, B)
    print(y, "\n")
    print("u x = y\nx = \n")
    x = GEPPLU(upper, y)
    print(x)

    end = time.time()
    timex = end - start
    print("time = ", timex)

    if (flag == 0):
        return 0 ,0
    elif (flag == 1):
        end=time.time()
        timex=end-start
        return x, timex
    else:
        return 2 ,2

###############################
## Gauss Elimination using partial pivoting
def GEPPLU(A, B):
    # output: x is the solution of Ax=b.
    print (A,"\n",B ,"\n\n")
    n =  len(A)
    if len(B) != n:
        raise ValueError("Invalid argument: incompatible sizes between A & B.", B.size, n)
    B = array(B)
    A = array(A)
    # parial pivoting
    for k in range(n-1):
        #Choose largest pivot element below
        maxindex = abs(A[k:,k]).argmax() + k
        if A[maxindex, k] == 0:
            raise ValueError("Matrix is singular.")
        #Swap rows
        if maxindex != k:
            A[[k,maxindex]] = A[[maxindex, k]]
            B[[k,maxindex]] = B[[maxindex, k]]
            print(A)
            print(B)
        for row in range(k+1, n):
            multiplier = A[row][k]/A[k][k]
            #the only one in this column since the rest are zero
            A[row][k] = multiplier
            for col in range(k + 1, n):
                A[row][col] = A[row][col] - multiplier*A[k][col]
            #Equation solution column
            B[row] = B[row] - multiplier*B[k]
            print(A)
            print(B)

    x = zeros(n)
    k = n-1
    x[k] = B[k]/A[k,k]
    while k >= 0:
        x[k] = (B[k] - dot(A[k,k+1:],x[k+1:]))/A[k,k]
        k = k-1
    return x
## Gauss Elimination using partial pivoting
## Gauss Elimination using partial pivoting
def GEPP(A, b, n):
    # output: x is the solution of Ax=b.
    start = time.time()
    A = array(A)
    b = array(b)
    aug = np.zeros((n,n+1))

    for i in range(n):
        for j in range(n+1):
            if (j != n):
                aug[i][j] = A[i][j]
            else:
                aug[i][j] = b[i]
    print("aug = \n",aug)

    # RANK OF MATRCIS
    rankA = np.linalg.matrix_rank(A)
    rankAug = np.linalg.matrix_rank(aug)
    print("rankA = ",rankA,"\nrankAug = ",rankAug,"\n")
    if(rankA == rankAug):
        if(rankA == n):
            print("The system has a Unique solution")
            flag = 1
        else:
            print("The system has a infinte number of solution")
            flag = 2
    else:
        print("The system doesnot has any solution")
        flag = 0

    n =  len(A)
    if len(b) != n:
        raise ValueError("Invalid argument: incompatible sizes between A & b.", b.size, n)
    A=array(A)
    # parial pivoting
    for k in range(n-1):
        #Choose largest pivot element below
        maxindex = abs(A[k:,k]).argmax() + k
        if A[maxindex, k] == 0:
            raise ValueError("Matrix is singular.")
        #Swap rows
        if maxindex != k:
            A[[k,maxindex]] = A[[maxindex, k]]
            b[[k,maxindex]] = b[[maxindex, k]]
            print(A,"\n")
            print(b,"\n\n")
        for row in range(k+1, n):
            multiplier = A[row][k]/A[k][k]
            #the only one in this column since the rest are zero
            A[row][k] = multiplier
            print(A,"\n")
            print(b,"\n\n")
            for col in range(k + 1, n):
                A[row][col] = A[row][col] - multiplier*A[k][col]
            #Equation solution column
            b[row] = b[row] - multiplier*b[k]
            print(A,"\n")
            print(b,"\n\n")

    x = zeros(n)
    k = n-1
    x[k] = b[k]/A[k,k]
    while k >= 0:
        x[k] = (b[k] - dot(A[k,k+1:],x[k+1:]))/A[k,k]
        k = k-1
    print("X = ",x)

    end = time.time()
    timex = end - start

    print("time = ",timex)

    if (flag == 0):
        return 0 , 0
    elif(flag == 1):
        return x , timex
    else:
        return 2 ,2

##########################################################33333

def gaussJordan(aug, a, n):
    start = time.time()
    # Making numpy array of n size and initializing
    # to zero for storing solution vector
    x = np.zeros(n)

    for i in range(n):
        for j in range(n):
            a[i][j] = aug[i][j]
    aug = array(aug)
    print("A = \n", a)

    # RANK OF MATRCIS
    rankA = np.linalg.matrix_rank(a)
    rankAug = np.linalg.matrix_rank(aug)
    print("rankA = ", rankA, "\nrankAug = ", rankAug, "\n")
    if (rankA == rankAug):
        if (rankA == n):
            print("The system has a Unique solution")
            flag = 1
        else:
            print("The system has a infinte number of solution")
            flag = 2
    else:
        print("The system doesnot has any solution")
        flag = 0
    # aug = [[0 for col in range(n+1)] for row in range(n)]

    # Applying Gauss Jordan Elimination
    for i in range(n):
        # Choose largest pivot element below
        maxindex = abs(aug[i:, i]).argmax() + i
        if aug[maxindex, i] == 0:
            raise ValueError("Matrix is singular.")
        # Swap rows
        if maxindex != i:
            aug[[i, maxindex]] = aug[[maxindex, i]]
            print("after partial pivoting\n", aug, "\n")

        for j in range(n):
            if (i == j):
                if (aug[i][i] != 1):
                    if (aug[i][i] != 0):
                        jordan = aug[i][i]
                        for k in range(n + 1):
                            aug[i][k] = float(aug[i][k] / jordan)

        print(aug, "\n")
        for j in range(n):
            if i != j:
                ratio = aug[j][i] / aug[i][i]

                for k in range(n + 1):
                    aug[j][k] = aug[j][k] - ratio * aug[i][k]
                print(aug, "\n")

    # Obtaining Solution

    for i in range(n):
        x[i] = aug[i][n]

    # Displaying solution
    print('\nRequired solution is: ')
    for i in range(n):
        print('X%d = %f' % (i, x[i]), end='\t')



    if (flag == 0):
        return 0 ,0
    elif (flag == 1):
        end = time.time()
        timex = end - start
        return x, timex
    else:
        return 2 ,2





###############################################################################################33



def evaluate(equation,variable):
    equation = re.sub('[a-z]',str(variable), equation)
    return eval(equation)


def bracketing (algorithm,equation,xl, xu, eps,  num_of_iterations):
    if(algorithm==1):
        return bracketing_algorithms(bisection,equation,xl, xu, eps, num_of_iterations)
    else:
        return bracketing_algorithms(false_position, equation, xl, xu, eps, num_of_iterations)

# Implementing Bracketing Method
def bracketing_algorithms(function,equation,xl, xu, eps, num_of_iterations):
   if evaluate(equation,xl) * evaluate(equation,xu) >= 0:
     raise BaseException("Wrong interval!")
   if xl > xu:
     xl, xu = xu, xl
   iterations = 0
   xr = function(xl,xu,equation)
   xr_old = 0
   start = time.time()
   error=abs(xr - xr_old) / abs(xr)
   start = datetime.datetime.now()
   while error >= eps / 100 and num_of_iterations != iterations:
         xr_old = xr

         if evaluate(equation,xr) == 0:
             end = datetime.datetime.now()
             return xr, iterations,(start-end).total_seconds()*1000,error
         elif evaluate(equation,xl) * evaluate(equation,xr) < 0:
           xu = xr
         else:
           xl = xr
         xr = function(xl,xu,equation)
         iterations += 1
         error = abs(xr - xr_old) / abs(xr)
         print("iteration  %d : " % (iterations),"xr: " ,xr ,"lower_x",xl,"Upper_x",xu)
   end = datetime.datetime.now()
   return xr,evaluate(equation,xr), iterations,(start-end).total_seconds()*1000, error




# Implementing Fixed Point Iteration Method
def fixedPointIteration(expression_g,x0, eps, num_of_iterations):

    xr_old=0
    xr=x0
    iterations = 0
    start = datetime.datetime.now()
    print(start)
    error = abs(xr - xr_old) / abs(xr)
    while error >= eps / 100 and num_of_iterations != iterations:
        xr_old = xr
        xr=evaluate(expression_g,xr_old)
        if(xr==xr_old) :
            iterations += 1
            end = datetime.datetime.now()

            return xr, iterations,(start-end).total_seconds()*1000,error
        iterations += 1
        error = abs(xr - xr_old) / abs(xr)
        print("iteration  %d :  xr=" % (iterations),xr)
    end = datetime.datetime.now()
    print((start-end).total_seconds()*1000)
    return xr, iterations,(start-end).total_seconds()*1000,error



# def check_convergence(equation):
#      equation = re.sub('[a-z]', 'x', equation)
#      symbols = {'x': Symbol('x', True)}
#      func = parse_expr(equation, my_symbols)
#      if(evaluate(diff(my_func, my_symbols['x']),1) <1) :
#          return true
#      else :
#          return false


def bisection(xl,xu,equation):
  return(xl + xu) / 2

def false_position(x0,x1,equation) :
  return x0 - (x1-x0) * evaluate(equation,x0)/ (evaluate(equation,x1) - evaluate(equation,x0) )


def gauss_Seidil(coef,ans,initial_guess,eps,iterations) :
    list=initial_guess
    num_of_iterations=0
    count=0
    error = [None] * len(ans)
    start = datetime.datetime.now()
    while num_of_iterations != iterations:
        list=initial_guess.copy()
        count =0
        for i in range(len(ans)):
            initial_guess[i]=ans[i]
            for j in range(len(ans)):
                if i != j:
                    initial_guess[i] -= (coef[i][j] * initial_guess[j])
            initial_guess[i] /= coef[i][i]
        num_of_iterations+=1
        for i in range(len(list)):
            error[i]=abs(initial_guess[i] - list[i]) / abs(initial_guess[i])
        for i in range(len(list)) :
            if(error[i] <= eps / 100 ) :
                count+=1
            if(count==len(list)) :
                end = datetime.datetime.now()
                return initial_guess ,num_of_iterations,int((start-end).total_seconds()*1000000),error
        print("iteration  %d :" % (num_of_iterations),initial_guess)
    end = datetime.datetime.now()
    return initial_guess, num_of_iterations,(start-end).total_seconds()*1000,error


def jacobi_iterations(coef,ans,initial_guess,eps,iterations) :
    list_old=[None]*len(ans)
    num_of_iterations=0
    count=0
    start = datetime.datetime.now()
    error=[None]*len(ans)
    while num_of_iterations != iterations:
        count =0
        for i in range(len(ans)):
            list_old[i]=ans[i]
            for j in range(len(ans)):
                if i != j:
                    list_old[i] -= (coef[i][j] * initial_guess[j])
            list_old[i] /= coef[i][i]
        num_of_iterations+=1
        for i in range(len(list_old)):
           error[i]= abs(list_old[i] - initial_guess[i]) / abs(list_old[i])
        for i in range(len(list_old)) :
            if( error[i] <= eps / 100 ) :
                count+=1
            if(count==len(list_old)) :
                end = datetime.datetime.now()
                return list_old ,num_of_iterations,(start-end).total_seconds()*1000,error
        initial_guess = list_old.copy()
        print("iteration  %d :" % (num_of_iterations),initial_guess)
    end = datetime.datetime.now()
    return initial_guess, num_of_iterations,(start-end).total_seconds()*1000,error

def plot_function(g_x,bracketing) :
    x = np.linspace(-6, 6, 50)

    # Plot degree 2 Taylor polynomial
    y2 = eval(g_x)
    plt.plot(x, y2, 'r-.', 'Degree 2')
    if(not bracketing) :
        y4 = x
        plt.plot(x, y4, 'g:', 'Degree 4')
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, color='k')
    plt.show()






################################################################################################



# Defining Function
def f(x, exp):
    exp = exp.lower()

    exp = exp.replace('^', '**')
    exp = exp.replace('log', 'math.log10')

    exp = exp.replace('sin', 'math.sin')
    exp = exp.replace('cos', 'math.cos')
    exp = exp.replace('tan', 'math.tan')
    exp = exp.replace('sinh', 'math.sinh')
    exp = exp.replace('cosh', 'math.cosh')
    exp = exp.replace('tanh', 'math.tanh')

    exp = exp.replace('pi', 'math.pi')
    exp = exp.replace('e', 'math.e')

    return eval(exp)



# Implementing Newton Raphson Method

def newtonRaphson(exp, x0, e, N):
    start=time.time()
    print('\n\n*** NEWTON RAPHSON METHOD IMPLEMENTATION ***')
    step = 1
    flag = 1
    condition = True
    x = smp.symbols('x')
    while condition:
      #  if dfdx() == 0.0:
       #     print('Divide by zero error!')
        #    break

        df=smp.diff(exp,x)

        x1 = x0 - f(x0,exp) / f(x0,str(df))
        print('Iteration-%d, x1 = %0.6f and f(x1) = %0.6f' % (step, x1, f(x1,exp)))
        x0 = x1
        step = step + 1

        if step > N:
            flag = 0
            break

        condition = abs(f(x1,exp)) > e

    if flag == 1:
        end = time.time()
        timex=end - start
        return x1 ,timex
    else:
        print('\nNot Convergent.')




# Note: You can combine above three section like this
#exp = input('Enter the function: ')
#x0 = float(input('Enter Guess: '))
#e = float(input('Tolerable Error: '))
#N = int(input('Maximum Step: '))

# Starting Newton Raphson Method
# newtonRaphson(exp, x0, e, N)


########################################################################################

# Defining Function

def fun(x, exp):
    exp = exp.lower()

    exp = exp.replace('^', '**')
    exp = exp.replace('log', 'math.log10')

    exp = exp.replace('sin', 'math.sin')
    exp = exp.replace('cos', 'math.cos')
    exp = exp.replace('tan', 'math.tan')
    exp = exp.replace('sinh', 'math.sinh')
    exp = exp.replace('cosh', 'math.cosh')
    exp = exp.replace('tanh', 'math.tanh')

    exp = exp.replace('pi', 'math.pi')
    exp = exp.replace('e', 'math.e')

    return eval(exp)


# Implementing Secant Method

def secant(exp, x0, x1, e, N):
    start = time.time()
    step = 1
    condition = True
    while condition:
        if (fun(x0, exp) == fun(x1, exp)):
            print('Divide by zero error!')
            break

        x2 = x0 - (x1 - x0) * fun(x0, exp) / (fun(x1, exp) - fun(x0, exp))
        print('Iteration-%d, x2 = %0.6f and f(x2) = %0.6f' % (step, x2, fun(x2, exp)))
        x0 = x1
        x1 = x2
        step = step + 1

        if step > N:
            print('Not Convergent!')
            break

        condition = abs(fun(x2, exp)) > e
    end = time.time()
    timex = end - start
    return x2, timex




