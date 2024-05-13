import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import libsbml

# Define the function to create activation and inhibition matrices
def create_matrices(filename):
    nw = pd.read_excel(filename)
    num_of_nodes = len(nw['Nodes'])
    mact = np.zeros((num_of_nodes, num_of_nodes))
    minh = np.zeros((num_of_nodes, num_of_nodes))
    node_names = nw['Nodes']

    for i in range(num_of_nodes):
        act_in_line = [node.strip() for node in str(nw['Activators'][i]).split(',')]
        inh_in_line = [node.strip() for node in str(nw['Inhibitors'][i]).split(',')]
        mact[i, :] = [1 if node in act_in_line else 0 for node in node_names]
        minh[i, :] = [1 if node in inh_in_line else 0 for node in node_names]

    return mact, minh, node_names, num_of_nodes

# Read network and extract data
mact, minh, node_names, num_of_nodes = create_matrices('CRT.xlsx')

# Initial values and simulation setup
xinit = np.random.rand(num_of_nodes)
clamped = np.zeros(num_of_nodes)
time_span = (0, 30)

# Define the ODE system function
def ODESysFun(t, X, NumOfNodes, Mact, Minh, Clamped):
    gamma = np.ones(NumOfNodes)  # Decay constant of each node
    h = 10  # Steepness of activation
    f = np.zeros(NumOfNodes)

    for i in range(NumOfNodes):
        Ract = Mact[i, :]
        Rinh = Minh[i, :]
        sum_alpha_X = np.dot(Ract, X)
        sum_beta_X = np.dot(Rinh, X)
        sum_alpha = np.sum(Ract)
        sum_beta = np.sum(Rinh)

        if np.any(Rinh == 0) and np.any(Ract):
            w = ((1 + sum_alpha) / sum_alpha) * (sum_alpha_X / (1 + sum_alpha_X))
        elif np.any(Ract == 0) and np.any(Rinh):
            w = 1 - ((1 + sum_beta) / sum_beta) * (sum_beta_X / (1 + sum_beta_X))
        elif np.any(Ract) and np.any(Rinh):
            w = (((1 + sum_alpha) / sum_alpha) * (sum_alpha_X / (1 + sum_alpha_X))) * (1 - ((1 + sum_beta) / sum_beta) * (sum_beta_X / (1 + sum_beta_X)))
        else:
            w = 0

        f[i] = (-np.exp(0.5 * h) + np.exp(-h * (w - 0.5))) / ((1 - np.exp(0.5 * h)) * (1 + np.exp(-h * (w - 0.5)))) - gamma[i] * X[i]
        
        if Clamped[i] == 1:
            f[i] = 0

    return f

# Solve the ODE system
sol = solve_ivp(ODESysFun, time_span, xinit, args=(num_of_nodes, mact, minh, clamped), method='RK45')

# Plotting and Visualization
fig, axs = plt.subplots(2, 3, figsize=(15, 10))
groups = {
    'Pro-Inflammatory': [0, 1, 2],  # Adjusted for Python's zero-based index
    'Anti-Inflammatory': [3],
    'Growth factors': [9],
    'ECM Destruction': [0, 1],
    'ECM Synthesis': [2],
    'Hypertrophy': [49]  # Adjusted index
}

for idx, (title, indices) in enumerate(groups.items()):
    ax = axs.flat[idx]
    ax.plot(sol.t, sol.y[indices, :].T, linewidth=2)
    ax.set(title=title, xlabel='Artificial Time', ylabel='Normalized Presence')
    ax.legend([node_names[i] for i in indices])
    plt.savefig('fig1.png')

# Steady States Bar Plot
plt.figure(figsize=(14, 8))  # Width, height in inches

# Create a bar plot
plt.bar(range(num_of_nodes), sol.y[:, -1])

# Set the x-ticks and their labels with increased font size
plt.gca().set_xticks(range(num_of_nodes))
plt.gca().set_xticklabels(node_names, rotation=45, fontsize=10)  # Adjust fontsize as needed

# Set the title with increased font size
plt.title('Current Steady State', fontsize=14)

# Optionally increase y-axis label size if needed
plt.ylabel('Value', fontsize=12)

# Save the figure
plt.tight_layout()  # This often helps with fitting labels that otherwise might be cut off
plt.savefig('fig.png')
