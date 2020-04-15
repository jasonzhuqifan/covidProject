from seir_main import SEIR 
import unittest
import numpy as np
from sklearn.metrics import r2_score

# This is the unit test file to test the function in SEIR class.
# Run: `python3 -m unittest test_sier_main.py` to run the tests.
class TestSeir(unittest.TestCase):

	def test_forward(self):
		self.seir = SEIR()
		S = 10000
		E = 10
		I = 100
		D = 50
		C = 100
		param = [0.1,0.2,0.3,0.4,0.5]
		iteration = 50
		result = self.seir.forward(S,E,I,D,C,param,iteration)
		self.assertEqual(len(result[0]), 5)
		self.assertEqual(len(result), iteration)

	def test_log_loss(self):
		self.seir = SEIR()
		y_pred = np.asarray([[1,2,3,4], [1,2,3,4],[2,2,3,4]])
		y_true = np.asarray([[1,2,3,4], [1,3,2,4],[2,3,4,4]])
		result = self.seir.logloss(y_true, y_pred)
		self.assertEqual(round(result,6), 0.6204060)

	def test_r2_score(self):
		self.seir = SEIR()
		y_pred = [[1,2,3,4], [1,2,3,4],[1,1]]
		y_true = [[1,2,3,4], [1,3,2,4], [2,1]]
		for i in range(len(y_pred)):
			result = self.seir.r2_score(y_true[i], y_pred[i])
			self.assertEqual(result, r2_score(y_true[i], y_pred[i]))

	def test_get_column_data(self):
		self.seir = SEIR()
		rows = [[1,2,3], [1,2,3],[2,2,3],[4,5,6],[66,66,44]]
		result = self.seir.get_column_data(rows)
		self.assertEqual(len(result[0]), len(rows))
		self.assertEqual(len(result), len(rows[0]))

	def test_score(self):
		self.seir = SEIR()
		self.seir.P = (0.1,0.2,0.3,0.4,0.5)
		S = 10000
		E = 10
		I = 100
		D = 50
		C = 100
		Y = [[1,2,3], [1,2,3],[66,66,44]]
		result = self.seir.score(S,E,I,D,C,Y)
		self.assertEqual(len(result), 2)
