# plot_legendre.py

import csv
import matplotlib.pyplot as plt
import os

# Read the CSV file that the C++ program created
x = []
polys = [[] for _ in range(6)]      # P0 to P5
ders  = [[] for _ in range(6)]      # P0' to P5'

with open('build/legendre.csv', newline='') as f:
    reader = csv.reader(f)
    next(reader)                    # skip header
    for row in reader:
        x.append(float(row[0]))
        for k in range(6):
            polys[k].append(float(row[1 + 2*k]))
            ders[k].append(float(row[2 + 2*k]))

# Plot
plt.figure(figsize=(12, 9))

# Polynomials
plt.subplot(2, 1, 1)
for k in range(6):
    plt.plot(x, polys[k], label=f'$P_{k}(x)$', linewidth=2.5)
plt.title('Legendre Polynomials $P_0(x)$ – $P_5(x)$', fontsize=16)
plt.xlabel('$x$')
plt.ylabel('$P_n(x)$')
plt.xlim(-1, 1)
plt.grid(True, alpha=0.3)
plt.legend()

# Derivatives
plt.subplot(2, 1, 2)
for k in range(6):
    plt.plot(x, ders[k], label=f"$P_{k}'(x)$", linewidth=2.5)
plt.title('Derivatives computed with Automatic Differentiation', fontsize=16)
plt.xlabel('$x$')
plt.ylabel("$dP_n/dx$")
plt.xlim(-1, 1)
plt.grid(True, alpha=0.3)
plt.legend()

# Save in multiple nice formats (all in the same folder as the CSV!)
output_png = os.path.join(os.path.dirname('build/legendre.csv'), "legendre_polynomials.png")

plt.savefig(output_png, dpi=300, bbox_inches='tight', facecolor='white')

print(f"Plot saved as:")
print(f"   {output_png}")

plt.tight_layout()
plt.show()