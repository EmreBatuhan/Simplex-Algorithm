
# This method reads the file in the given path and builds up the starting variables we need
def parse_input_txt(input_path):
    input_file = open(input_path, "r")
    num_of_constraints, num_of_vars = map(int, input_file.readline().split())
    # All of the variables are constructed as a two dimensional array
    # Because we want all of them to be in matrix form
    c = [list(map(int, input_file.readline().split()))]
    A = []
    b = []
    for i in range(num_of_constraints): 
        line = list(map(int, input_file.readline().split()))
        A.append(line[:-1])
        b.append([line[-1]])
    return num_of_constraints, num_of_vars, A, b, c

# Since the problem is not given in canonical form we add slacks to each constraint
# Every constaint has the lesser than ot equal to constraint and the right-hand sides are non-negative
# So we simply add a slack to each constraint and the canonical form is achieved
def add_slacks(A, c, num_of_constraints): 
    # New A matrix is constructed by adding the corresponding 1s and 0s to the end
    for row in range(num_of_constraints):
        for i in range(num_of_constraints): 
            if i == row: 
                A[row].append(1)
            else: 
                A[row].append(0)
    # Similarly the objective function variables -c- should be extended with 0s
    for i in range(num_of_constraints):
        c[0].append(0)
    return A, c, len(A[0])

def matrix_multiply(A, B): 
    row_A = len(A)
    col_A = len(A[0])
    col_B = len(B[0])
    result = [[0 for i in range(col_B)] for j in range(row_A)]
    for i in range(row_A): 
        for j in range(col_B): 
            res = 0
            for k in range(col_A): 
                res += A[i][k] * B[k][j]
            result[i][j] = res
    return result

def transposeMatrix(A):
    return list(map(list,zip(*A)))

def getMatrixMinor(A,i,j):
    return [row[:j] + row[j+1:] for row in (A[:i]+A[i+1:])]

def getMatrixDeterminant(A):
    #base case for 2x2 matrix
    if len(A) == 2:
        return A[0][0]*A[1][1]-A[0][1]*A[1][0]

    determinant = 0
    for c in range(len(A)):
        determinant += ((-1)**c)*A[0][c]*getMatrixDeterminant(getMatrixMinor(A,0,c))
    return determinant

#This method finds the inverse of a given matrix using determinant and the adjoint matrix of A
def inverse(A):
    determinant = getMatrixDeterminant(A)
    #special case for 2x2 matrix:
    if len(A) == 2:
        return [[A[1][1]/determinant, -1*A[0][1]/determinant],
                [-1*A[1][0]/determinant, A[0][0]/determinant]]

    #find matrix of cofactors
    cofactors = []
    for r in range(len(A)):
        cofactorRow = []
        for c in range(len(A)):
            minor = getMatrixMinor(A,r,c)
            cofactorRow.append(((-1)**(r+c)) * getMatrixDeterminant(minor))
        cofactors.append(cofactorRow)
    cofactors = transposeMatrix(cofactors)
    for r in range(len(cofactors)):
        for c in range(len(cofactors)):
            cofactors[r][c] = cofactors[r][c]/determinant
    return cofactors

def matrix_sub(A, B): 
    row = len(A)
    col = len(A[0])
    result = [[0 for i in range(col)] for j in range(row)]
    for i in range(row): 
        for j in range(col):
            result[i][j] = A[i][j] - B[i][j]
    return result

#This method takes a list of base indexes that should be included in the B matrix and A itself
def build_base_matrix(A, base_list): 
    result = []
    for row in A: 
        result_row = []
        for i in base_list: 
            result_row.append(row[i])
        result.append(result_row)
    return result

#This method checks the final objective function variables and returns the the index of the smallest value
#If there are no variable with a negative value this method returns -1 which corresponds to not found
def control_objective_row(C): 
    min_el = 0
    min_ind = -1
    for ind in range(len(C[0])):
        el = C[0][ind]
        if el < min_el: 
            min_el = el
            min_ind = ind
    return min_ind

#In order to advance in the Simplex algorithm an incoming index should be replaced with a chosen index
def find_outgoing_index(A, b, incoming_index):
    search_col = incoming_index
    min_ratio = float("inf")
    min_ind = -1
    #A ratio test is done to determine an outgoing index
    for i in range(len(A)):
        #If the variable inside the chosen column is not positive it cant be chosen
        #The small value instead of a direct 0 is due to float presicion
        if A[i][search_col] <= 0.0000001:
            continue
        ratio = b[i][0] / A[i][search_col]
        #Searching for the smallest non-negative ratio
        if ratio >= 0 and ratio < min_ratio:
            min_ind = i
            min_ratio = ratio
    #There is a possibility that this method returns -1
    #Which means that there are no constraints limiting this incoming variable
    #Therefore the problem is unbounded
    return min_ind

#Prints out the current instance of the tableu
def print_simplex_tableu(A, b, c, z): 
    for i in range(len(A)):
        for j in range(len(A[0])):
            print("{:6.2f}".format(A[i][j]),end=" ")
        print("| {:6.2f}".format(b[i][0]))
    print("-"*7*(len(A[0])+2))
    for i in range(len(c[0])):
        print("{:6.2f}".format(c[0][i]),end=" ")
    print("| {:6.2f}".format(z[0][0]))
    print()

input_path = "Data3.txt"
num_of_constraints, num_of_vars, A, b, c = parse_input_txt(input_path)
A, c, num_of_vars = add_slacks(A, c, num_of_constraints)
z0 = [[0]]

#The starting base is consisting of the last elements due to the slack variables we added
base = [num_of_vars + i - num_of_constraints for i in range(num_of_constraints)]

#The star variables are the loop variables
#They are given the initial values here and they will be changed repeatedly in the loop
A_star = A
c_star = c
b_star = b
z0_star = z0
while True:
    #Printing out the current tableu
    print_simplex_tableu(A_star, b_star, c_star, z0_star)

    #Trying to improve current tableu by adding a variable to the base
    base_incoming_index = control_objective_row(c_star)
    if base_incoming_index == -1:
        #There are no variables in the objective function which choosing it improves current situation
        #Current tableu is optimal
        initial_num_of_vars = num_of_vars-num_of_constraints # The count of initial variables that are not slacks
        solution = [0 for i in range(initial_num_of_vars)]
        for i in range(len(base)):
            if base[i] < initial_num_of_vars:
                solution[base[i]] = b_star[i][0]

        print("Solution : [",end="")
        for i in range(len(solution)):
            print("{:.2f} ".format(solution[i]),end="")
        print("]")
        print("Optimal Result: {:.2f}".format(z0_star[0][0]*-1))
        break
    #Since we add a variable to the base, we should remove one variable from the base using ratio test
    base_outgoing_index = find_outgoing_index(A_star, b_star, base_incoming_index)
    if base_outgoing_index == -1:
        #Ratio test could not find any feasible point
        #Which means the problem is unbounded
        print("NO SOLUTION")
        break

    #The swapping part of incoming and outgoing variables
    base[base_outgoing_index]=base_incoming_index

    #Since the base is changed, the A,b,c and z0 matrices are updated
    B = build_base_matrix(A_star, base)
    c_B = build_base_matrix(c_star, base)
    B_inverse = inverse(B)

    A_star = matrix_multiply(B_inverse, A_star)
    b_star = matrix_multiply(B_inverse, b_star)
    c_star = matrix_sub(c_star, matrix_multiply(c_B, A_star))
    z0_star = matrix_sub(z0_star, matrix_multiply(c_B, b_star))

