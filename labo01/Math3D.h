#pragma once

/**
 * @file Math3D.h
 *
 * @brief Fonctions pour l'intinialisation et la manipulation de matrices de
 * rotation, des matrices de transformations en coordonnées homogènes et les
 * vecteurs 3D.
 *
 * Nom: Andre Godlove Hessede
 * Code permanent : HESA89050004
 * Email : andre-godlove.hessede.1@ens.etsmtl.ca
 *
 */

#include "Matrix.h"
#include "Vector.h"
#include "Operators.h"

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <math.h>
#include <iostream>


namespace gti320 {

	// Deux types de vecteurs 3D considérés ici
	typedef Vector<double, 3> Vector3d;
	typedef Vector<float, 3> Vector3f;

	// Dans le cadre de ce projet, nous considérons seulement deux
	// cas :
	//
	//  - les rotations
	//  - les translations
	//
	// Deux types de matrices en coordonnées homogèes :
	typedef Matrix<double, 4, 4, ColumnStorage> Matrix4d;
	typedef Matrix<float, 4, 4, ColumnStorage> Matrix4f;
	// 
	// Deux types de matrices pour les rotations
	typedef Matrix<double, 3, 3, ColumnStorage> Matrix3d;
	typedef Matrix<float, 3, 3, ColumnStorage> Matrix3f;

	/**
	 * Initialise et retourne la matrice identité
	 */
	template<>
	inline void Matrix4d::setIdentity()
	{
		// TODO affecter la valeur 0.0 partout, sauf sur la diagonale principale où c'est 1.0.
		//      Note: ceci est une redéfinition d'une fonction membre!
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				(*this)(i, j) = 0.0;
			}
		}

		for (int i = 0; i < 4; ++i) {
			(*this)(i, i) = 1.0;
		}
	}

	/**
	 * Calcul de la matrice inverse, SPÉCIALISÉ pour le cas d'une matrice de
	 * transformation en coordonnées homogènes.
	 *
	 * TODO (vous pouvez supposer qu'il s'agit d'une matrice de transformation
	 * en coordonnées homogènes)
	 */
	template<>
	inline Matrix4d Matrix4d::inverse() const
	{
		// TODO : implémenter
		Matrix4d result;
		result.setZero();

		Matrix3d rotation;
		Vector3d translation;
		Matrix3d rotationInverse;
		Vector3d translationInverse;

		rotation = this->block(0, 0, 3, 3);
		translation = this->block(0, 3, 3, 1);
		rotationInverse = rotation.transpose<double, 3, 3, ColumnStorage>();
		translationInverse = (-1.0 * translation);
		translation = rotationInverse * translationInverse;

		result.block(0, 0, 3, 3) = rotationInverse;
	
		for (int i = 0; i < translation.size(); ++i) {
			result(i, 3) = translation(i);
		}

		result(3, 3) = 1.0;


		return result; 
	}

	/**
	 * Calcul de la matrice inverse, SPÉCIALISÉ pour le cas d'une matrice de rotation.
	 *
	 * (vous pouvez supposer qu'il s'agit d'une matrice de rotation)
	 */
	template<>
	inline Matrix3d Matrix3d::inverse() const
	{
		return this->transpose<double, 3, 3, ColumnStorage>();
	}


	/**
	 * Multiplication d'une matrice 4x4 avec un vecteur 3D où la matrice
	 * représente une transformation en coordonnées homogène.
	 */
	template <typename _Scalar>
	Vector<_Scalar, 3> operator*(const Matrix<_Scalar, 4, 4, ColumnStorage>& A, const Vector<_Scalar, 3>& v)
	{
		// TODO : implémenter

		Vector<_Scalar, 3> result;
		result.setZero();

		// Calcul de la transformation homogène
		for (int i = 0; i < 3; ++i) {
			result(i) = A(i, 0) * v(0) + A(i, 1) * v(1) + A(i, 2) * v(2) + A(i, 3);
		}

		return result;
	}

	/**
	 * 
	 * Création de la matrice 3d de mise à echelle
	 */
	template <typename _Scalar>
	Matrix3d scaling(const _Scalar scale)
	{
		// Implémentation supplémentaire

		Matrix3d result;
		result.setZero();

		for (int i = 0; i < 3; ++i) {
			result(i, i) = scale;
		}

		return result;
	}

	/**
	 *
	 * Produit vectoriel de deux vecteur 3d
	 */
	template <typename _Scalar>
	Vector3d cross(const Vector<_Scalar, 3>& v1, const Vector<_Scalar, 3>& v2)
	{
		// Implémentation supplémentaire
		Vector<_Scalar, 3> result;
		result.setZero();

		result(0) = v1(1) * v2(2) - v1(2) * v2(1);
		result(1) = v1(2) * v2(0) - v1(0) * v2(2);
		result(2) = v1(0) * v2(1) - v1(1) * v2(0);

		return result;
	}


	/**
	 * Initialise et retourne la matrice de rotation définie par les angles
	 * d'Euler XYZ exprimés en radians.
	 *
	 * La matrice doit correspondre au produit : Rz*Ry*Rx.
	 */
	template<typename _Scalar>
	static Matrix<_Scalar, 3, 3> makeRotation(_Scalar x, _Scalar y, _Scalar z)
	{
		// TODO : implémenter
		_Scalar cx = std::cos(x);
		_Scalar sx = std::sin(x);
		_Scalar cy = std::cos(y);
		_Scalar sy = std::sin(y);
		_Scalar cz = std::cos(z);
		_Scalar sz = std::sin(z);

		Matrix<_Scalar, 3, 3> rotation;
		rotation.setZero();

		// Matrice de rotation autour de l'axe X (Roll)
		rotation(0, 0) = 1.0;
		rotation(0, 1) = 0.0;
		rotation(0, 2) = 0.0;
		rotation(1, 0) = 0.0;
		rotation(1, 1) = cx;
		rotation(1, 2) = -sx;
		rotation(2, 0) = 0.0;
		rotation(2, 1) = sx;
		rotation(2, 2) = cx;

		// Matrice de rotation autour de l'axe Y (Pitch)
		Matrix<_Scalar, 3, 3> rotationY;
		rotationY(0, 0) = cy;
		rotationY(0, 1) = 0.0;
		rotationY(0, 2) = sy;
		rotationY(1, 0) = 0.0;
		rotationY(1, 1) = 1.0;
		rotationY(1, 2) = 0.0;
		rotationY(2, 0) = -sy;
		rotationY(2, 1) = 0.0;
		rotationY(2, 2) = cy;

		// Matrice de rotation autour de l'axe Z (Yaw)
		Matrix<_Scalar, 3, 3> rotationZ;
		rotationZ(0, 0) = cz;
		rotationZ(0, 1) = -sz;
		rotationZ(0, 2) = 0.0;
		rotationZ(1, 0) = sz;
		rotationZ(1, 1) = cz;
		rotationZ(1, 2) = 0.0;
		rotationZ(2, 0) = 0.0;
		rotationZ(2, 1) = 0.0;
		rotationZ(2, 2) = 1.0;

		// Produit des trois matrices de rotation
		rotation = rotationZ * (rotationY * rotation);

		return rotation;
	}

}
