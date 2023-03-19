#!/usr/bin/python3
import math

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1


class GeneSequencing:

	def __init__(self):
		pass

	# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
	# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
	# how many base pairs to use in computing the alignment

	def align(self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length

		if banded:

			# d and k designate our bandwidth – d is the number of cells to calculate from the diagonal,
			# 	k is the total bandwidth
			d = 3
			k = 7

			table_2d = []
			prev_2d = []

			POINT_LEFT = 0
			POINT_TOP = 1
			POINT_TOPLEFT = 2
			POINT_TOPRIGHT = 3

			# This initializes the first row of the table
			table_2d.append([])
			prev_2d.append([])
			counter = 0

			# This initializes the first banded row
			for i in range(d + 1):

				if i > d:
					break
				table_2d[0].append(counter)
				counter += 5
				prev_2d[0].append(POINT_LEFT)

			# Set to 5 since it shares 0 with the first row
			counter = 5

			# i is our rows, so it has to cover the entire length of the string
			for i in range(align_length + 1):

				# This check covers the first two sequences (polynomial and exponential)
				if len(seq2) - len(seq1) > 1000:
					too_big = True
					break
				else:
					too_big = False

				if i > len(seq2) or i > len(seq1) + 1:
					break

				if i == 0:
					continue

				table_2d.append([])
				prev_2d.append([])

				# We use this counter to access the correct chars in seq1 (the string that goes down the rows)
				# 	We need it because we have to access 3 chars less than and 3 chars greater than the diagonal
				# 	So it will run from -3 to 3 and be added to i inside the second for loop (j)
				character_counter = - 3

				# This is how we switch the number of columns that we need at the beginning and end of the matrix
				if i >= align_length - d:
					columns_left = align_length - i + d + 1
				elif i >= len(seq1) - d:
					columns_left = len(seq1) - i + d + 1
				else:
					columns_left = math.inf

				if i <= d:
					columns_to_check = d + i + 1
				else:
					columns_to_check = math.inf

				# This loop keeps our time and space to O(kn) instead of making it O(mn)
				for j in range(k):

					# Break when there's no more columns –– Helps manage row length
					if columns_to_check == 0:
						break
					columns_to_check -= 1

					# This initializes the first column of the table
					if j == 0 and i <= d:
						table_2d[i].append(counter)
						counter += 5
						prev_2d[i].append(POINT_TOP)
						continue

					# This if is for the first 4 rows, the else for all the others.
					# 	The way I manage this, we don't change whether diagonals are called top_left or not.
					# 	They just get treated differently in the prev array
					if i <= d:

						if columns_to_check > 1:
							top = table_2d[i - 1][j] + 5
						else:
							top = math.inf
						left = table_2d[i][j - 1] + 5
						top_left = table_2d[i - 1][j - 1]

						if seq1[j - 1] == seq2[i - 1]:
							top_left -= 3
						else:
							top_left += 1

					else:

						if j == 0:
							left = math.inf
						else:
							left = table_2d[i][j - 1] + 5
						if j == k - 1:
							top = math.inf
						else:
							top = table_2d[i - 1][j + 1] + 5
						top_left = table_2d[i - 1][j]

						# The character_counter is how we know where in seq1 to access, since we only iterate through
						# 	k times
						if seq1[i + character_counter - 1] == seq2[i - 1]:
							top_left -= 3
						else:
							top_left += 1

					# priority –– left, top, diagonal –– we need the i <= d checks to know which way to point
					if left <= top:
						if left <= top_left:
							table_2d[i].append(left)
							prev_2d[i].append(POINT_LEFT)
						else:
							table_2d[i].append(top_left)
							if i <= d:
								prev_2d[i].append(POINT_TOPLEFT)
							else:
								prev_2d[i].append(POINT_TOP)
					elif top <= top_left:
						table_2d[i].append(top)
						if i <= d:
							prev_2d[i].append(POINT_TOP)
						else:
							prev_2d[i].append(POINT_TOPRIGHT)
					else:
						table_2d[i].append(top_left)
						if i <= d:
							prev_2d[i].append(POINT_TOPLEFT)
						else:
							prev_2d[i].append(POINT_TOP)

					# this should never be more than 3 –– it's how we know when to stop iterating if there's >7 columns
					if i >= align_length - d or i >= len(seq2) - d:
						columns_left -= 1
						if columns_left == 0:
							break

					character_counter += 1
			# This is how we handle big discrepencies in sequence length
			if too_big:
				score = math.inf
				alignment1 = "No alignment possible"
				alignment2 = "No alignment possible"

			else:
				alignment1 = ""
				alignment2 = ""
				# seq1 counter allows us to access the string since our column counter stays between
				# 	0 < column_counter < 7. row_counter also functions as seq2 counter
				row_counter = len(table_2d) - 2
				column_counter = len(table_2d[-1]) - 1
				if len(seq1) > align_length:
					seq1_counter = align_length - 1
				else:
					seq1_counter = len(seq1) - 1

				# Start with the initial character in the bottom right corner
				alignment1 = seq1[seq1_counter] + alignment1
				alignment2 = seq2[row_counter] + alignment2

				while True:

					if row_counter == 0 and column_counter == 0:
						break

					if row_counter <= d:
						if prev_2d[row_counter][column_counter] == POINT_LEFT:
							column_counter -= 1
							seq1_counter -= 1
							alignment1 = seq1[seq1_counter] + alignment1
							alignment2 = '-' + alignment2
						elif prev_2d[row_counter][column_counter] == POINT_TOP:
							row_counter -= 1
							alignment1 = '-' + alignment1
							alignment2 = seq2[row_counter] + alignment2
						else:
							column_counter -= 1
							seq1_counter -= 1
							row_counter -= 1
							alignment1 = seq1[seq1_counter] + alignment1
							alignment2 = seq2[row_counter] + alignment2
					else:
						if prev_2d[row_counter][column_counter] == POINT_LEFT:
							column_counter -= 1
							seq1_counter -= 1
							alignment1 = seq1[seq1_counter] + alignment1
							alignment2 = '-' + alignment2
						elif prev_2d[row_counter][column_counter] == POINT_TOPRIGHT:
							row_counter -= 1
							column_counter += 1
							alignment1 = '-' + alignment1
							alignment2 = seq2[row_counter] + alignment2
						else: # this functions as the diagonal sub/match now
							seq1_counter -= 1
							row_counter -= 1
							alignment1 = seq1[seq1_counter] + alignment1
							alignment2 = seq2[row_counter] + alignment2

				score = table_2d[-1][-1]
				alignment1 = alignment1[:100]
				alignment2 = alignment2[:100]

		else:
			# initialize table – first value is rows, second value is columns
			table_2d = []
			prev_2d = []

			# Values for the prev array:
			POINT_LEFT = 0
			POINT_TOP = 1
			POINT_TOPLEFT = 2

			# This initializes the first row of the table
			table_2d.append([])
			prev_2d.append([])
			counter = 0
			for i in range(align_length + 1):
				if i > len(seq1):
					break
				table_2d[0].append(counter)
				counter += 5
				prev_2d[0].append(POINT_LEFT)

			counter = 5
			for i in range(align_length + 1):

				if i > len(seq2):
					break

				if i == 0:
					continue

				table_2d.append([])
				prev_2d.append([])

				for j in range(align_length + 1):

					if j > len(seq1):
						break

					# This initializes the first column of the table
					if j == 0:
						table_2d[i].append(counter)
						counter += 5
						prev_2d[i].append(POINT_TOP)
						continue

					# comparing values within the table as they get populated
					left = table_2d[i][j - 1] + 5
					top = table_2d[i - 1][j] + 5
					top_left = table_2d[i - 1][j - 1]
					if seq1[j - 1] == seq2[i - 1]:
						top_left -= 3
					else:
						top_left += 1

					# Priority order: left, top, top-left
					if left <= top:
						if left <= top_left:
							table_2d[i].append(left)
							prev_2d[i].append(POINT_LEFT)
						else:
							table_2d[i].append(top_left)
							prev_2d[i].append(POINT_TOPLEFT)
					elif top <= top_left:
						table_2d[i].append(top)
						prev_2d[i].append(POINT_TOP)
					else:
						table_2d[i].append(top_left)
						prev_2d[i].append(POINT_TOPLEFT)

			alignment1 = ""
			alignment2 = ""
			row_counter = len(table_2d) - 2
			column_counter = len(table_2d[0]) - 2

			# Start with the initial character in the bottom right corner
			alignment1 = seq1[column_counter] + alignment1
			alignment2 = seq2[row_counter] + alignment2

			# Traverse the prev array from end back to beginning to build the strings
			while True:

				if row_counter == 0 and column_counter == 0:
					break

				if prev_2d[row_counter][column_counter] == POINT_LEFT:
					column_counter -= 1
					alignment1 = seq1[column_counter] + alignment1
					alignment2 = '-' + alignment2
				elif prev_2d[row_counter][column_counter] == POINT_TOP:
					row_counter -= 1
					alignment1 = '-' + alignment1
					alignment2 = seq2[row_counter] + alignment2
				else:
					column_counter -= 1
					row_counter -= 1
					alignment1 = seq1[column_counter] + alignment1
					alignment2 = seq2[row_counter] + alignment2

			score = table_2d[-1][-1]
			alignment1 = alignment1[:100]
			alignment2 = alignment2[:100]

		return {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}
