#pragma once

/**
 * @file Operators.h
 *
 * @brief Opérateurs arithmétiques pour les matrices et les vecteurs.
 *
 * Nom: Andre Godlove Hessede
 * Code permanent : HESA89050004
 * Email : andre-godlove.hessede.1@ens.etsmtl.ca
 *
 */

#include "Matrix.h"
#include "Vector.h"

 /**
  * Implémentation de divers opérateurs arithmétiques pour les matrices et les vecteurs.
  */
namespace gti320 {

	/**
	 * Multiplication : Matrice * Matrice (générique)
	 */
	template <typename _Scalar, int RowsA, int ColsA, int StorageA, int RowsB, int ColsB, int StorageB>
	Matrix<_Scalar, RowsA, ColsB> operator*(const Matrix<_Scalar, RowsA, ColsA, StorageA>& A, const Matrix<_Scalar, RowsB, ColsB, StorageB>& B)
	{
		// TODO implémenter
		assert(A.cols() == B.rows());

		Matrix<_Scalar, RowsA, ColsB> result(A.rows(), B.cols());
		result.setZero();

		for (int i = 0; i < A.rows(); ++i) {
			for (int j = 0; j < B.cols(); ++j) {
				_Scalar sum = 0;
				for (int k = 0; k < A.cols(); ++k) {
					sum += A(i, k) * B(k, j);
				}
				result(i, j) = sum;
			}
		}
		return result;
	}

	/**
	 * Multiplication : Matrice (colonne) * Matrice (ligne)
	 *
	 * Spécialisation de l'opérateur de multiplication pour le cas où les matrices
	 * ont un stockage à taille dynamique et où la matrice de gauche utilise un
	 * stockage par colonnes et celle de droite un stockage par lignes.
	 */
	template <typename _Scalar>
	Matrix<_Scalar, Dynamic, Dynamic> operator*(const Matrix<_Scalar, Dynamic, Dynamic, ColumnStorage>& A, const Matrix<_Scalar, Dynamic, Dynamic, RowStorage>& B)
	{
		// TODO : implémenter
		assert(A.cols() == B.rows());

		Matrix<_Scalar, Dynamic, Dynamic> result(A.rows(), B.cols());
		result.setZero();

		for (int k = 0; k < A.cols(); ++k) {
			for (int j = 0; j < B.cols(); ++j) {
				for (int i = 0; i < A.rows(); ++i) {
					result(i, j) += A(i, k) * B(k, j);
				}
			}
		}

		return result;
	}

	/**
	 * Multiplication : Matrice (ligne) * Matrice (colonne)
	 *
	 * Spécialisation de l'opérateur de multiplication pour le cas où les matrices
	 * ont un stockage à taille dynamique et où la matrice de gauche utilise un
	 * stockage par lignes et celle de droite un stockage par colonnes.
	 */
	template <typename _Scalar>
	Matrix<_Scalar, Dynamic, Dynamic> operator*(const Matrix<_Scalar, Dynamic, Dynamic, RowStorage>& A, const Matrix<_Scalar, Dynamic, Dynamic, ColumnStorage>& B)
	{
		// TODO : implémenter
		assert(A.cols() == B.rows());

		Matrix<_Scalar, Dynamic, Dynamic> result(A.rows(), B.cols());
		result.setZero();

		for (int i = 0; i < A.rows(); ++i) {
			for (int j = 0; j < B.cols(); ++j) {
				for (int k = 0; k < A.cols(); ++k) {
					result(i, j) += A(i, k) * B(k, j);
				}
			}
		}

		return result;
	}


	/**
	 * Addition : Matrice + Matrice (générique)
	 */
	template <typename _Scalar, int Rows, int Cols, int StorageA, int StorageB>
	Matrix<_Scalar, Rows, Cols> operator+(const Matrix<_Scalar, Rows, Cols, StorageA>& A, const Matrix<_Scalar, Rows, Cols, StorageB>& B)
	{
		// TODO : implémenter
		Matrix<_Scalar, Rows, Cols> result(A.rows(), B.cols());
		result.setZero();

		for (int i = 0; i < A.rows(); ++i) {
			for (int j = 0; j < B.cols(); ++j) {
				result(i, j) = A(i, j) + B(i, j);
			}
		}

		return result;
	}

	/**
	 * Addition : Matrice (colonne) + Matrice (colonne)
	 *
	 * Spécialisation de l'opérateur d'addition pour le cas où les deux matrices
	 * sont stockées par colonnes.
	 */
	template <typename _Scalar>
	Matrix<_Scalar, Dynamic, Dynamic> operator+(const Matrix<_Scalar, Dynamic, Dynamic, ColumnStorage>& A, const Matrix<_Scalar, Dynamic, Dynamic, ColumnStorage>& B)
	{
		// TODO : implémenter
		assert(A.rows() == B.rows() && A.cols() == B.cols());

		int numRows = A.rows();
		int numCols = A.cols();

		Matrix<_Scalar, Dynamic, Dynamic> result(numRows, numCols);
		result.setZero();

		for (int j = 0; j < numCols; ++j) {
			for (int i = 0; i < numRows; ++i) {
				result(i, j) = A(i, j) + B(i, j);
			}
		}
		return result;
	}

	/**
	 * Addition : Matrice (ligne) + Matrice (ligne)
	 *
	 * Spécialisation de l'opérateur d'addition pour le cas où les deux matrices
	 * sont stockées par lignes.
	 */
	template <typename _Scalar>
	Matrix<_Scalar, Dynamic, Dynamic, RowStorage> operator+(const Matrix<_Scalar, Dynamic, Dynamic, RowStorage>& A, const Matrix<_Scalar, Dynamic, Dynamic, RowStorage>& B)
	{
		// TODO : implémenter
		assert(A.rows() == B.rows() && A.cols() == B.cols());

		int numRows = A.rows();
		int numCols = A.cols();

		Matrix<_Scalar, Dynamic, Dynamic, RowStorage> result(numRows, numCols);
		result.setZero();

		for (int i = 0; i < numRows; ++i) {
			for (int j = 0; j < numCols; ++j) {
				result(i, j) = A(i, j) + B(i, j);
			}
		}

		return result;
	}

	/**
	 * Multiplication  : Scalaire * Matrice (colonne)
	 *
	 * Spécialisation de l'opérateur de multiplication par un scalaire pour le
	 * cas d'une matrice stockée par colonnes.
	 */
	template <typename _Scalar, int _Rows, int _Cols>
	Matrix<_Scalar, _Rows, _Cols, ColumnStorage> operator*(const _Scalar& a, const Matrix<_Scalar, _Rows, _Cols, ColumnStorage>& A)
	{
		// TODO : implémenter
		Matrix<_Scalar, _Rows, _Cols, ColumnStorage> result(A.rows(), A.cols());
		result.setZero();

		for (int j = 0; j < A.cols(); ++j) {
			for (int i = 0; i < A.rows(); ++i) {
				result(i, j) = a * A(i, j);
			}
		}
		return result;
	}

	/**
	 * Multiplication  : Scalaire * Matrice (ligne)
	 *
	 * Spécialisation de l'opérateur de multiplication par un scalaire pour le
	 * cas d'une matrice stockée par lignes.
	 */
	template <typename _Scalar, int _Rows, int _Cols>
	Matrix<_Scalar, _Rows, _Cols, RowStorage> operator*(const _Scalar& a, const Matrix<_Scalar, _Rows, _Cols, RowStorage>& A)
	{
		// TODO : implémenter
		Matrix<_Scalar, _Rows, _Cols, RowStorage> result(A.rows(), A.cols());
		result.setZero();

		for (int i = 0; i < A.rows(); ++i) {
			for (int j = 0; j < A.cols(); ++j) {
				result(i, j) = a * A(i, j);
			}
		}
		return result;
	}

	/**
	 * Multiplication : Matrice (ligne) * Vecteur
	 *
	 * Spécialisation de l'opérateur de multiplication matrice*vecteur pour le
	 * cas où la matrice est représentée par lignes.
	 */
	template <typename _Scalar, int _Rows, int _Cols>
	Vector<_Scalar, _Rows> operator*(const Matrix<_Scalar, _Rows, _Cols, RowStorage>& A, const Vector<_Scalar, _Cols>& v)
	{
		// TODO : implémenter
		assert(A.cols() == v.size());

		Vector<_Scalar, _Rows> result(A.cols());
		result.setZero();

		for (int i = 0; i < A.rows(); ++i) {
			_Scalar sum = 0;
			for (int j = 0; j < A.cols(); ++j) {
				sum += A(i, j) * v(j);
			}
			result(i) = sum;
		}

		return result;
	}

	/**
	 * Multiplication : Matrice (colonne) * Vecteur
	 *
	 * Spécialisation de l'opérateur de multiplication matrice*vecteur pour le
	 * cas où la matrice est représentée par colonnes.
	 */
	template <typename _Scalar, int _Rows, int _Cols>
	Vector<_Scalar, _Rows> operator*(const Matrix<_Scalar, _Rows, _Cols, ColumnStorage>& A, const Vector<_Scalar, _Cols>& v)
	{
		// TODO : implémenter
		assert(A.cols() == v.size());

		Vector<_Scalar, _Rows> result(A.cols());
		result.setZero();

		for (int j = 0; j < A.cols(); ++j) {
			for (int i = 0; i < A.rows(); ++i) {
				result(i) += A(i, j) * v(j);
			}
		}

		return result;
	}

	/**
	 * Multiplication : Scalaire * Vecteur
	 */
	template <typename _Scalar, int _Rows>
	Vector<_Scalar, _Rows> operator*(const _Scalar& a, const Vector<_Scalar, _Rows>& v)
	{
		// TODO : implémenter
		Vector<_Scalar, _Rows> result(v.size());
		result.setZero();

		for (int i = 0; i < v.size(); ++i) {
			result(i) = a * v(i);
		}

		return result;
	}


	/**
	 * Addition : Vecteur + Vecteur
	 */
	template <typename _Scalar, int _RowsA, int _RowsB>
	Vector<_Scalar, _RowsA> operator+(const Vector<_Scalar, _RowsA>& a, const Vector<_Scalar, _RowsB>& b)
	{
		// TODO : implémenter
		assert(a.size() == b.size());

		Vector<_Scalar, _RowsA> result(a.size());
		result.setZero();

		for (int i = 0; i < a.size(); ++i) {
			result(i) = a(i) + b(i);
		}

		return result;
	}

	/**
	 * Soustraction : Vecteur - Vecteur
	 */
	template <typename _Scalar, int _RowsA, int _RowsB>
	Vector<_Scalar, _RowsA> operator-(const Vector<_Scalar, _RowsA>& a, const Vector<_Scalar, _RowsB>& b)
	{
		// TODO : implémenter
		assert(a.size() == b.size());

		Vector<_Scalar, _RowsA> result(a.size());
		result.setZero();

		for (int i = 0; i < a.size(); ++i) {
			result(i) = a(i) - b(i);
		}

		return result;
	}
}
