import pandas as pd
import matplotlib.pyplot as plt


stats = pd.read_csv(
	'results/statistics.csv',
	names=[ 'iteration',
			'n_living_cells',
			'average_cell_energy',
			'average_cell_chem',
			'average_cell_toxin',
			'total_environmental_chem',
			'total_environmental_toxin'
	]
)

stats.iloc[:, 0:5].plot(x='iteration')
plt.xlabel('Iteration number')
plt.savefig('results/CellStatistics.png', dpi=300)

stats.iloc[:, [0, 5, 6]].plot(x='iteration')
plt.xlabel('Iteration number')
plt.savefig('results/EnvironmentStatistics.png', dpi=300)