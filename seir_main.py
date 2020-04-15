import numpy as np
from scipy.optimize import dual_annealing, minimize
from collections import namedtuple

SEIR_PARAM = namedtuple('SEIRparm', ['alpha1', 'alpha2', 'beta', 'sigma', 'gamma'])

class SEIR(object):

	def __init__(self, P=None):

		self.P = P
	# Each iteration to find the next set of parameters with modified SEIR model.
	def forward(self, S, E, I, D, C, param, max_iter):
		a1, a2, b, s, g = param
		est = []
		for t in range(max_iter):
			S_ = S - a1 * E - a2 * I + s * I
			E_ = E + a1 * E + a2 * I - b * E
			I_ = I + b * E - s * I - g * I
			D_ = D + g * I
			C_ = C + s * I
			S, E, I, D, C = S_, E_, I_, D_, C_
			est.append([S, E, I, D, C])
		return est

	# Metric function that computes the log loss of real data and predicted data.
	def logloss(self, obs, est):
		assert len(obs) == len(est)
		loss = ((np.log2(obs + 1) - np.log2(est + 1)) ** 2).sum()
		return loss

	# Metric function that computes the r2 score of real data and predicted data.
	def r2_score(self, y_true, y_pred):
		y_true = np.asarray(y_true)
		y_pred = np.asarray(y_pred)
		numerator = ((y_true - y_pred) ** 2).sum(axis=0,
														  dtype=np.float64)
		denominator = ((y_true - np.average(
			y_true, axis=0, weights=None)) ** 2).sum(axis=0,
															  dtype=np.float64)
		return 1 - numerator / denominator

	# Helper function that returns a list of data as a form of [confirmed, death, recovered].
	# Input: data  -  a list of three lists containing confirmed data, death data and recovered data in each list.
	def get_IDC(self, data):
		res = []
		for i in range(len(data)):
			res.append([data[i][2], data[i][3], data[i][4]])
		return res

	# Optimized function that uses log loss to fit the model
	def opt_fn(self, param, s, e, i, d, c, obs):
		est = self.forward(s, e, i, d, c, param, len(obs))
		return self.logloss(obs, np.asarray(self.get_IDC(est)))

	# Fit the input data and find the best parameters with optimized function defined above
	def fit(self, initS, initE, initI, initD, initC, Y):
		args = (initS, initE, initI, initD, initC, np.asarray(Y))
		param = [(0, 1), ] * 5
		result = dual_annealing(self.opt_fn, param, args=args, seed=30, maxiter=10)['x']
		self.P = SEIR_PARAM(*result)

	# Convert rows data to three-column data (confirmed, death, recovered).
	def get_column_data(self, rows):
		col1 = []
		col2 = []
		col3 = []
		for i in range(len(rows)):
			col1.append(rows[i][0])
			col2.append(rows[i][1])
			col3.append(rows[i][2])
		return col1, col2, col3

	# Computer log loss and average r2 score.
	def score(self, init_S, init_E, init_I, init_D, init_C, Y):
		est = self.get_IDC(self.predict(init_S, init_E, init_I, init_D, init_C, len(Y)))
		loss = self.logloss(np.asarray(Y), np.asarray(est))

		y1, y2, y3 = self.get_column_data(Y)
		est1, est2, est3 = self.get_column_data(est)
		r1 = self.r2_score(y1, est1)
		r2 = self.r2_score(y2, est2)
		r3 = self.r2_score(y3, est3)

		return loss, (r1 + r2 + r3) / 3

	# Make predictions for T days with initial values.
	def predict(self, init_S, init_E, init_I, init_D, init_C, T):
		return self.forward(init_S, init_E, init_I, init_D, init_C, self.P, T)

# Convert training data as column form.
def get_train_data(confirmed_data, death_data, recovered_data):
	train = []
	for i in range(len(confirmed_data)):
		train.append([confirmed_data[i], death_data[i], recovered_data[i]])
	return train

# Main function that finds the best params (alpha1, alpha2, beta, sigma, gamma) and optimal initial E.
def search_param(init_S, init_E, step, confirmed, death, recovered):
	seir = SEIR()
	train = get_train_data(confirmed, death, recovered)
	min_loss, max_r2, best_param, likeli_potential = float('inf'), float('-inf'), None, 0
	for potential in range(0, int(init_E), int(step)):
		seir.fit(int(init_S), potential, train[0][0], train[0][1], train[0][2], train)
		loss, r2 = seir.score(int(init_S), potential, train[0][0], train[0][1],  train[0][2], Y=train)
		if loss < min_loss and r2 > max_r2:
			min_loss, max_r2, best_param, likeli_potential = loss, r2, seir.P, potential
	seir.P = best_param
	return best_param, likeli_potential
