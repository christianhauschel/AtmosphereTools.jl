using DifferentialEquations

# function solve_tf(N, D, u0, t_end)

u0 = [0.0]
N = [1.0, 1.0]
D = [1.0]
t_end = 1.0

# Order of the numerator and denominator
n = length(N) - 1
m = length(D) - 1




# Build system n first-order ODEs
# function ode(dy, y, p, t)
#     dy[1] = y[2]
#     dy[2] = y[3]
#     ...
#     dy[m-1] = y[m]
#     dy[m] = 1 / D[m] * (
#         N[n] * 
#     )
# end 


# prob = ODEProblem(ode, u0, t_end)
# println(prob)
# sol = solve(prob, u0)
# return sol
# return nothing
# end


# sol = solve_tf(N, D, u0, t_end)