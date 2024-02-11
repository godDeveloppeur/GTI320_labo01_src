#pragma once

/**
 * @file Vector.h
 *
 * @brief Implémentation de vecteurs simples
 *
 * Nom: Andre Godlove Hessede
 * Code permanent : HESA89050004
 * Email : andre-godlove.hessede.1@ens.etsmtl.ca
 *
 */

#include <cmath>
#include "MatrixBase.h"

namespace gti320
{

	/**
	 * Classe vecteur générique.
	 *
	 * Cette classe réutilise la classe `MatrixBase` et ses spécialisations de
	 * templates pour les manipulation bas niveau.
	 */
	template <typename _Scalar = double, int _Rows = Dynamic>
	class Vector : public MatrixBase<_Scalar, _Rows, 1>
	{
	public:
		/**
		 * Constructeur par défaut
		 */
		Vector() : MatrixBase<_Scalar, _Rows, 1>() {}

		/**
		 * Contructeur à partir d'un taille (rows).
		 */
		explicit Vector(int rows) : MatrixBase<_Scalar, _Rows, 1>(rows, 1) {}

		/**
		 * Constructeur de copie
		 */
		Vector(const Vector& other) : MatrixBase<_Scalar, _Rows, 1>(other) {}

		/**
		 * Destructeur
		 */
		~Vector() {}

		/**
		 * Opérateur de copie
		 */
		Vector& operator=(const Vector& other)
		{
			// TODO implémenter
			MatrixBase<_Scalar, _Rows, 1>::operator=(other);
			return *this;
		}

		/**
		 * Opérateur de copie à partir d'une sous-matrice
		 */
		template <typename _OtherScalar, int OtherRows, int _OtherCols, int _OtherStorage>
		Vector& operator=(const SubMatrix<_OtherScalar, OtherRows, _OtherCols, _OtherStorage>& submatrix)
		{
			// Implémentation supplémentaire
			assert(submatrix.cols() == 1);
			for (int i = 0; i < MatrixBase<_Scalar, _Rows, 1>::size(); ++i)
			{
				MatrixBase<_Scalar, _Rows, 1>::m_storage.data()[i] = submatrix(i, 0);
			}
			return *this;
		}


		/**
		 * Opérateur de copie à partir d'une liste
		 */
		Vector& operator=(const _Scalar* arr)
		{
			// Implémentation supplémentaire
			for (int i = 0; i < MatrixBase<_Scalar, _Rows, 1>::size(); i++) {
				MatrixBase<_Scalar, _Rows, 1>::m_storage.data()[i] = arr[i];
			}
			return *this;
		}

		/**
		 * Accesseur à une entrée du vecteur (lecture seule)
		 */
		_Scalar operator()(int i) const
		{
			// TODO implémenter
			return MatrixBase<_Scalar, _Rows, 1>::m_storage.data()[i];
		}

		/**
		 * Accesseur à une entrée du vecteur (lecture et écriture)
		 */
		_Scalar& operator()(int i)
		{
			// TODO implémenter
			return MatrixBase<_Scalar, _Rows, 1>::m_storage.data()[i];
		}

		/**
		 * Modifie le nombre de lignes du vecteur
		 */
		void resize(int _rows)
		{
			MatrixBase<_Scalar, _Rows, 1>::resize(_rows, 1);
		}

		/**
		 * Produit scalaire de *this et other.
		 */
		inline _Scalar dot(const Vector& other) const
		{
			// TODO implémenter
			assert(this->size() == other.size());
			_Scalar result = 0.0;
			for (int i = 0; i < other.size(); ++i)
			{
				result += MatrixBase<_Scalar, _Rows, 1>::m_storage.data()[i] * other(i);
			}
			return result;
		}

		/**
		 * Retourne la norme euclidienne du vecteur
		 */
		inline _Scalar norm() const
		{
			// TODO implémenter
			_Scalar result = 0.0;
			for (int i = 0; i < MatrixBase<_Scalar, _Rows, 1>::size(); ++i)
			{
				result += pow(MatrixBase<_Scalar, _Rows, 1>::m_storage.data()[i], 2);
			}

			return sqrt(result);
		}

		/**
		 * Retourne la norme de manhattan du vecteur
		 */
		inline _Scalar normManhattan() const
		{
			// Implémentation supplémentaire
			_Scalar result = 0.0;
			for (int i = 0; i < MatrixBase<_Scalar, _Rows, 1>::size(); ++i)
			{
				result += abs(MatrixBase<_Scalar, _Rows, 1>::m_storage.data()[i]);
			}

			return result;
		}

		/**
		 * Retourne la norme infinie du vecteur
		 */
		inline _Scalar normInfinie() const
		{
			// Implémentation supplémentaire
			_Scalar result = 0.0;
			for (int i = 0; i < MatrixBase<_Scalar, _Rows, 1>::size(); ++i)
			{
				result = std::max(result, abs(MatrixBase<_Scalar, _Rows, 1>::m_storage.data()[i]));
			}

			return result;
		}

		/**
		 * Retourne la norme P du vecteur
		 */
		inline _Scalar normP(const int& p) const
		{
			// Implémentation supplémentaire
			_Scalar result = 0.0;
			for (int i = 0; i < MatrixBase<_Scalar, _Rows, 1>::size(); ++i)
			{
				result += std::pow(abs(MatrixBase<_Scalar, _Rows, 1>::m_storage.data()[i]), p);
			}

			return std::pow(result, 1.0 / p);
		}
	};
}
