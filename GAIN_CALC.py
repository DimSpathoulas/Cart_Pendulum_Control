import numpy as np
from scipy.linalg import solve_continuous_are

# system parameters
M = 1.0     # Cart mass
m = 0.1     # Pendulum mass
l = 0.5     # Pendulum length
g = 9.81    # Gravitational acceleration

# state-space matrices A and B (linearized system)
A = np.array([[0, 0, 1, 0],
              [0, 0, 0, 1],
              [0, -m*g/M, 0, 0],
              [0, (M+m)*g/(M*l), 0, 0]])

B = np.array([[0],
              [0],
              [1/M],
              [-1/(M*l)]])

# cost matrices Q and R
Q = np.diag([10.0, 10.0, 1.0, 1.0])  # State cost (penalizing state errors)
R = np.array([[0.1]])  # Control cost (penalizing control effort)

# continuous-time algebraic Riccati equation (ARE)
P = solve_continuous_are(A, B, Q, R)

# optimal feedback gain matrix K
K = np.linalg.inv(R) @ B.T @ P

print("The computed LQR gain K is:")
print(K)
