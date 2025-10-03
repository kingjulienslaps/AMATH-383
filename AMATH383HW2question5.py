import numpy as np
import matplotlib.pyplot as plt

# Define parameters
mu = np.linspace(0.5, 2, 500)  # Range of mu values
theta_branch_stable = lambda mu: np.sqrt(6 * (mu - 1) / mu)  # Stable branches for theta
theta_branch_unstable = lambda mu: -np.sqrt(6 * (mu - 1) / mu)  # Unstable branches

# Mask values where branches don't exist (mu < 1)
mu_stable = mu[mu > 1]
theta_stable = theta_branch_stable(mu_stable)
theta_unstable = theta_branch_unstable(mu_stable)

# Create the plot
plt.figure(figsize=(8, 6))
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Theta = 0 line
plt.axvline(1, color='red', linestyle='--', label=r'$\mu = 1$ (bifurcation point)')  # Bifurcation point

# Plot stable and unstable branches
plt.plot(mu_stable, theta_stable, label='Stable Branch', color='blue', linewidth=2)
plt.plot(mu_stable, theta_unstable, color='blue', linewidth=2)
plt.scatter(1, 0, color='red', zorder=5, label='Bifurcation Point')

# Labeling the diagram
plt.title('Bifurcation Diagram', fontsize=14)
plt.xlabel(r'$\mu$', fontsize=12)
plt.ylabel(r'$\theta$', fontsize=12)
plt.legend(fontsize=10)
plt.grid(alpha=0.4)
plt.ylim(-2, 2)
plt.xlim(0.5, 2)
plt.tight_layout()
plt.show()
