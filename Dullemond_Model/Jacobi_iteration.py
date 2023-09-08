import numpy as np

# Example from the video:
# from youtube https://www.youtube.com/watch?v=zPDp_ewoyhM
def jacobian_example(xy):
    x, y = xy
    return [[1, 2],
            [2*x, 8*y]]

def function_example(xy):
    x, y = xy
    return [x + 2*y - 2, x**2 + 4*y**2 - 4]

# From the exercise:
def function_exercise(xyz):
    x, y, z = xyz
    return [x + y + z - 3,
            x**2 + y**2 + z**2 - 5,
            np.exp(x) + x*y - x*z - 1]

def jacobian_exercise(xyz):
    x, y, z = xyz
    return [[1, 1, 1],
            [2*x, 2*y, 2*z],
            [np.exp(x) + y - z, x, -x]]



def iterative_newton(fun, x_init, jacobian):
    max_iter = 50
    epsilon = 1e-8

    x_last = x_init

    for k in range(max_iter):
        # Solve J(xn)*( xn+1 - xn ) = -F(xn):
        J = np.array(jacobian(x_last))
        F = np.array(fun(x_last))

        diff = np.linalg.solve( J, -F )
        x_last = x_last + diff

        # Stop condition:
        if np.linalg.norm(diff) < epsilon:
            print('convergence!, nre iter:', k )
            break

    else: # only if the for loop end 'naturally'
        print('not converged')

    return x_last

# For the exercice:
x_sol = iterative_newton(function_exercise, [2.0,1.0,2.0], jacobian_exercise)
print('solution exercice:', x_sol )
print('F(sol)', function_exercise(x_sol) )

# For the example:
x_sol = iterative_newton(function_example, [1.0,2.0], jacobian_example)
print('solution example:', x_sol )
print( function_example(x_sol) )