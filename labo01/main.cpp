/**
 * @file main.cpp
 *
 * @brief Unit tests for a simple linear algebra library.
 *
 * Nom: Andre Godlove Hessede
 * Code permanent : HESA89050004
 * Email : andre-godlove.hessede.1@ens.etsmtl.ca
 *
 */

#include "Matrix.h"
#include "Vector.h"
#include "Math3D.h"
#include "Operators.h"

#include <gtest/gtest.h>
#include <chrono>

using namespace gti320;

/**
 * Multiplication  matrice * vecteur,  utilisant une implémentation naive
 */
template<typename _Scalar>
static inline Vector<_Scalar, Dynamic> naiveMatrixMult(const Matrix<_Scalar, Dynamic, Dynamic, ColumnStorage>& A, const Vector<_Scalar, Dynamic>& v)
{
	assert(A.cols() == v.rows());

	Vector<_Scalar, Dynamic> b(A.rows());
	assert(b.rows() == A.rows());

	for (int i = 0; i < A.rows(); ++i) {
		b(i) = 0.0;
		for (int j = 0; j < A.cols(); ++j) {
			b(i) += A(i, j) * v(j);
		}
	}

	return b;
}

/**
 * Addition  matrice + matrice,  utilisant une implémentation naive
 */
template<typename _Scalar>
static inline Matrix<_Scalar, Dynamic, Dynamic, ColumnStorage> naiveMatrixAddition(const Matrix<_Scalar, Dynamic, Dynamic, ColumnStorage>& A, const Matrix<_Scalar, Dynamic, Dynamic, ColumnStorage>& B)
{
	assert(A.cols() == B.cols() && A.rows() == B.rows());

	Matrix<_Scalar, Dynamic, Dynamic, ColumnStorage> C(A.rows(), A.cols());
	assert(C.rows() == A.rows() && C.cols() == A.cols());
	for (int i = 0; i < C.rows(); ++i) {
		for (int j = 0; j < C.cols(); ++j) {
			C(i, j) = A(i, j) + B(i, j);
		}
	}
	return C;
}

/**
 * Multiplication  matrice * matrice,  utilisant une implémentation naive.
 */
template<typename _Scalar, int _Storage>
static inline Matrix<_Scalar, Dynamic, Dynamic, _Storage> naiveMatrixMult(const Matrix<_Scalar, Dynamic, Dynamic, _Storage>& A, const Matrix<_Scalar, Dynamic, Dynamic, _Storage>& B)
{
	assert(A.cols() == B.rows());
	Matrix<_Scalar, Dynamic, Dynamic> product(A.rows(), B.cols());
	for (int i = 0; i < A.rows(); ++i)
	{
		for (int j = 0; j < B.cols(); ++j)
		{
			for (int k = 0; k < A.cols(); ++k)
			{
				product(i, j) += A(i, k) * B(k, j);
			}
		}
	}
	return product;
}

// Test les matrice avec redimensionnement dynamique
TEST(TestLabo1, DynamicMatrixTests)
{
	// Crée une matrice à taille dynamique
	// (note : les valeurs par défaut du patron de la classe `Matrix` mettent le
	// le nombre de ligne et de colonnes à `Dynamic`)
	Matrix<double> M(3, 5);
	EXPECT_EQ(M.cols(), 5);
	EXPECT_EQ(M.rows(), 3);

	// Redimensionne la matrice
	M.resize(100, 1000);
	EXPECT_EQ(M.cols(), 1000);
	EXPECT_EQ(M.rows(), 100);

	// Test - stockage par colonnes
	Matrix<double, Dynamic, Dynamic, ColumnStorage> ColM(100, 100);
	ColM.setZero();
	ColM(0, 0) = 1.0;
	ColM(99, 99) = 99.0;
	ColM(10, 33) = 5.0;
	EXPECT_EQ(ColM(0, 0), 1.0);
	EXPECT_EQ(ColM(10, 33), 5.0);
	EXPECT_EQ(ColM(99, 99), 99.0);

	// Test - stockage par lignes
	Matrix<double, Dynamic, Dynamic, RowStorage> RowM(5, 4);
	RowM.setZero();
	RowM(0, 0) = 2.1;
	RowM(3, 3) = -0.2;
	RowM(4, 3) = 1.2;
	EXPECT_EQ(RowM.rows(), 5);
	EXPECT_EQ(RowM.cols(), 4);
	EXPECT_DOUBLE_EQ(RowM(0, 0), 2.1);
	EXPECT_DOUBLE_EQ(RowM(3, 3), -0.2);
	EXPECT_DOUBLE_EQ(RowM(4, 3), 1.2);
	EXPECT_DOUBLE_EQ(RowM(3, 2), 0.0);

	// Transposée
	const auto RowMT = RowM.transpose();
	EXPECT_EQ(RowMT.rows(), 4);
	EXPECT_EQ(RowMT.cols(), 5);
	EXPECT_DOUBLE_EQ(RowMT(0, 0), 2.1);
	EXPECT_DOUBLE_EQ(RowMT(3, 3), -0.2);
	EXPECT_DOUBLE_EQ(RowMT(3, 4), 1.2);
	EXPECT_DOUBLE_EQ(RowMT(2, 3), 0.0);
}



/**
 * Test pour les vecteurs à taille dynamique
 */
TEST(TestLabo1, DynamicVectorSizeTest)
{
	Vector<double> v(5);
	v.setZero();

	EXPECT_EQ(v.rows(), 5);

	v.resize(3);
	EXPECT_EQ(v.rows(), 3);

	v(0) = 1.0;
	v(1) = 2.0;
	v(2) = 3.0;

	EXPECT_DOUBLE_EQ(v.norm(), 3.7416573867739413855837487323165);

	Vector<double, Dynamic> v2(3);
	v2.setZero();
	v2(1) = 2.0;

	EXPECT_DOUBLE_EQ(v2.dot(v), 4.0);
	EXPECT_DOUBLE_EQ(v2(0), 0.0);
	EXPECT_DOUBLE_EQ(v2(1), 2.0);
	EXPECT_DOUBLE_EQ(v2(2), 0.0);
}

/**
 * Test pour les matrice à taille fixe
 */
TEST(TestLabo1, Matrix4x4SizeTest)
{
	Matrix4d M;
	M.setZero();

	EXPECT_EQ(M.cols(), 4);
	EXPECT_EQ(M.rows(), 4);
}

/**
 * Test pour les opérateurs d'arithmétique matricielle.
 */
TEST(TestLabo1, MatrixMatrixOperators)
{
	// Opérations arithmétiques avec matrices à taille dynamique
	{
		// Test : matrice identité
		Matrix<double> A(6, 6);
		A.setIdentity();
		EXPECT_DOUBLE_EQ(A(0, 0), 1.0);
		EXPECT_DOUBLE_EQ(A(1, 1), 1.0);
		EXPECT_DOUBLE_EQ(A(2, 2), 1.0);
		EXPECT_DOUBLE_EQ(A(3, 3), 1.0);
		EXPECT_DOUBLE_EQ(A(4, 4), 1.0);
		EXPECT_DOUBLE_EQ(A(5, 5), 1.0);
		EXPECT_DOUBLE_EQ(A(0, 1), 0.0);
		EXPECT_DOUBLE_EQ(A(1, 0), 0.0);

		// Test : produit  scalaire * matrice
		const double alpha = 2.5;
		Matrix<double> B = alpha * A;
		EXPECT_DOUBLE_EQ(B(0, 0), alpha);
		EXPECT_DOUBLE_EQ(B(1, 1), alpha);
		EXPECT_DOUBLE_EQ(B(2, 2), alpha);
		EXPECT_DOUBLE_EQ(B(3, 3), alpha);
		EXPECT_DOUBLE_EQ(B(4, 4), alpha);
		EXPECT_DOUBLE_EQ(B(5, 5), alpha);
		EXPECT_DOUBLE_EQ(B(0, 1), 0.0);
		EXPECT_DOUBLE_EQ(B(1, 0), 0.0);

		// Test : produit  matrice * matrice
		Matrix<double> C = A * B;
		EXPECT_DOUBLE_EQ(C(0, 0), A(0, 0) * B(0, 0));
		EXPECT_DOUBLE_EQ(C(1, 1), A(1, 1) * B(1, 1));
		EXPECT_DOUBLE_EQ(C(2, 2), A(2, 2) * B(2, 2));
		EXPECT_DOUBLE_EQ(C(3, 3), A(3, 3) * B(3, 3));
		EXPECT_DOUBLE_EQ(C(4, 4), A(4, 4) * B(4, 4));
		EXPECT_DOUBLE_EQ(C(5, 5), A(5, 5) * B(5, 5));
		EXPECT_DOUBLE_EQ(C(0, 1), 0.0);
		EXPECT_DOUBLE_EQ(C(2, 3), 0.0);

		// Test : addition  matrice * matrice
		Matrix<double> A_plus_B = A + B;
		EXPECT_DOUBLE_EQ(A_plus_B(0, 0), A(0, 0) + B(0, 0));
		EXPECT_DOUBLE_EQ(A_plus_B(1, 1), A(1, 1) + B(1, 1));
		EXPECT_DOUBLE_EQ(A_plus_B(2, 2), A(2, 2) + B(2, 2));
		EXPECT_DOUBLE_EQ(A_plus_B(3, 3), A(3, 3) + B(3, 3));
		EXPECT_DOUBLE_EQ(A_plus_B(4, 4), A(4, 4) + B(4, 4));
		EXPECT_DOUBLE_EQ(A_plus_B(5, 5), A(5, 5) + B(5, 5));
		EXPECT_DOUBLE_EQ(A_plus_B(0, 1), 0.0);
		EXPECT_DOUBLE_EQ(A_plus_B(2, 3), 0.0);
	}

	// Opérations arithmétique avec matrices à stockage par lignes et par
	// colonnes.
	{
		// Création d'un matrice à stockage par lignes
		Matrix<double, Dynamic, Dynamic, RowStorage> A(5, 5);
		A(0, 0) = 0.8147;    A(0, 1) = 0.0975;    A(0, 2) = 0.1576;    A(0, 3) = 0.1419;    A(0, 4) = 0.6557;
		A(1, 0) = 0.9058;    A(1, 1) = 0.2785;    A(1, 2) = 0.9706;    A(1, 3) = 0.4218;    A(1, 4) = 0.0357;
		A(2, 0) = 0.1270;    A(2, 1) = 0.5469;    A(2, 2) = 0.9572;    A(2, 3) = 0.9157;    A(2, 4) = 0.8491;
		A(3, 0) = 0.9134;    A(3, 1) = 0.9575;    A(3, 2) = 0.4854;    A(3, 3) = 0.7922;    A(3, 4) = 0.9340;
		A(4, 0) = 0.6324;    A(4, 1) = 0.9649;    A(4, 2) = 0.8003;    A(4, 3) = 0.9595;    A(4, 4) = 0.6787;

		// Test : transposée (le résultat est une matrice à stockage par
		//        colonnes)
		Matrix<double, Dynamic, Dynamic, ColumnStorage> B = A.transpose();

		// Test : multiplication  matrix(ligne) * matrice(colonne)
		// Note : teste seulement la première et la dernière colonne
		const auto C = A * B;
		EXPECT_NEAR(C(0, 0), 1.14815820000000, 1e-3); EXPECT_NEAR(C(0, 4), 1.31659795000000, 1e-3);
		EXPECT_NEAR(C(1, 0), 1.00133748000000, 1e-3); EXPECT_NEAR(C(1, 4), 2.04727044000000, 1e-3);
		EXPECT_NEAR(C(2, 0), 0.99433707000000, 1e-3); EXPECT_NEAR(C(2, 4), 2.82896409000000, 1e-3);
		EXPECT_NEAR(C(3, 0), 1.63883925000000, 1e-3); EXPECT_NEAR(C(3, 4), 3.28401323000000, 1e-3);
		EXPECT_NEAR(C(4, 0), 1.31659795000000, 1e-3); EXPECT_NEAR(C(4, 4), 3.35271580000000, 1e-3);


		// Test : multiplication  matrice(colonne) * matrice(ligne)
		// Note : teste seulement la première et la dernière colonne
		const auto C2 = B * A;
		EXPECT_NEAR(C2(0, 0), 2.73456805000000, 1e-3); EXPECT_NEAR(C2(0, 4), 1.95669703000000, 1e-3);
		EXPECT_NEAR(C2(1, 0), 1.88593811000000, 1e-3); EXPECT_NEAR(C2(1, 4), 2.08742862000000, 1e-3);
		EXPECT_NEAR(C2(2, 0), 2.07860468000000, 1e-3); EXPECT_NEAR(C2(2, 4), 1.94727447000000, 1e-3);
		EXPECT_NEAR(C2(3, 0), 1.94434955000000, 1e-3); EXPECT_NEAR(C2(3, 4), 2.27675041000000, 1e-3);
		EXPECT_NEAR(C2(4, 0), 1.95669703000000, 1e-3); EXPECT_NEAR(C2(4, 4), 2.48517748000000, 1e-3);

		// Test : addition  matrice(ligne) + matrice(ligne)
		// Note : teste seulement la première et la dernière colonne
		const auto A_plus_A = A + A;
		EXPECT_DOUBLE_EQ(A_plus_A(0, 0), A(0, 0) + A(0, 0)); EXPECT_DOUBLE_EQ(A_plus_A(0, 4), A(0, 4) + A(0, 4));
		EXPECT_DOUBLE_EQ(A_plus_A(1, 0), A(1, 0) + A(1, 0)); EXPECT_DOUBLE_EQ(A_plus_A(1, 4), A(1, 4) + A(1, 4));
		EXPECT_DOUBLE_EQ(A_plus_A(2, 0), A(2, 0) + A(2, 0)); EXPECT_DOUBLE_EQ(A_plus_A(2, 4), A(2, 4) + A(2, 4));
		EXPECT_DOUBLE_EQ(A_plus_A(3, 0), A(3, 0) + A(3, 0)); EXPECT_DOUBLE_EQ(A_plus_A(3, 4), A(3, 4) + A(3, 4));
		EXPECT_DOUBLE_EQ(A_plus_A(4, 0), A(4, 0) + A(4, 0)); EXPECT_DOUBLE_EQ(A_plus_A(4, 4), A(4, 4) + A(4, 4));

		// Test : addition  matrice(colonne) + matrice(colonne)
		// Note : teste seulement la première et la dernière colonne
		const auto B_plus_B = B + B;
		EXPECT_DOUBLE_EQ(B_plus_B(0, 0), B(0, 0) + B(0, 0)); EXPECT_DOUBLE_EQ(B_plus_B(0, 4), B(0, 4) + B(0, 4));
		EXPECT_DOUBLE_EQ(B_plus_B(1, 0), B(1, 0) + B(1, 0)); EXPECT_DOUBLE_EQ(B_plus_B(1, 4), B(1, 4) + B(1, 4));
		EXPECT_DOUBLE_EQ(B_plus_B(2, 0), B(2, 0) + B(2, 0)); EXPECT_DOUBLE_EQ(B_plus_B(2, 4), B(2, 4) + B(2, 4));
		EXPECT_DOUBLE_EQ(B_plus_B(3, 0), B(3, 0) + B(3, 0)); EXPECT_DOUBLE_EQ(B_plus_B(3, 4), B(3, 4) + B(3, 4));
		EXPECT_DOUBLE_EQ(B_plus_B(4, 0), B(4, 0) + B(4, 0)); EXPECT_DOUBLE_EQ(B_plus_B(4, 4), B(4, 4) + B(4, 4));

	}
}

/**
 * Test pour la multiplication  matrice * vecteur
 */

TEST(TestLabo1, MatrixVectorOperators)
{
	// Vecteur à taille dynamique
	Vector<double> v(5);
	v(0) = 1.0;
	v(1) = 2.0;
	v(2) = 4.0;
	v(3) = 8.0;
	v(4) = 16.0;

	// Test : multiplication par la matrice identité
	{
		Matrix<double> M(5, 5);
		M.setIdentity();

		const auto b = M * v;
		EXPECT_DOUBLE_EQ(b(0), 1.0);
		EXPECT_DOUBLE_EQ(b(1), 2.0);
		EXPECT_DOUBLE_EQ(b(2), 4.0);
		EXPECT_DOUBLE_EQ(b(3), 8.0);
		EXPECT_DOUBLE_EQ(b(4), 16.0);
	}

	// Test : multiplication par une matrice à taille dynamique avec stockage par ligne.
	{
		Matrix<double, Dynamic, Dynamic, RowStorage> M(5, 5);
		M.setIdentity();
		M = 2.0 * M;

		Vector<double> b2 = M * v;
		EXPECT_DOUBLE_EQ(b2(0), 2.0);
		EXPECT_DOUBLE_EQ(b2(1), 4.0);
		EXPECT_DOUBLE_EQ(b2(2), 8.0);
		EXPECT_DOUBLE_EQ(b2(3), 16.0);
		EXPECT_DOUBLE_EQ(b2(4), 32.0);
	}
}

/**
 * Opérateurs d'arithmétique vectorielle
 */


TEST(TestLabo1, VectorOperators)
{
	Vector<double> v(5);
	v(0) = 0.1;
	v(1) = 0.2;
	v(2) = 0.4;
	v(3) = 0.8;
	v(4) = 1.6;

	// Test : multiplication  scalaire * vecteur
	const double alpha = 4.0;
	const auto v2 = alpha * v;
	EXPECT_DOUBLE_EQ(v2(0), alpha * v(0));
	EXPECT_DOUBLE_EQ(v2(1), alpha * v(1));
	EXPECT_DOUBLE_EQ(v2(2), alpha * v(2));
	EXPECT_DOUBLE_EQ(v2(3), alpha * v(3));
	EXPECT_DOUBLE_EQ(v2(4), alpha * v(4));

	// Test : addition  vecteur + vecteur
	const auto v3 = v + v2;
	EXPECT_DOUBLE_EQ(v3(0), v(0) + v2(0));
	EXPECT_DOUBLE_EQ(v3(1), v(1) + v2(1));
	EXPECT_DOUBLE_EQ(v3(2), v(2) + v2(2));
	EXPECT_DOUBLE_EQ(v3(3), v(3) + v2(3));
	EXPECT_DOUBLE_EQ(v3(4), v(4) + v2(4));
}


/**
 * Mathématiques 3D
 */

TEST(TestLabo1, Math3D)
{
	// Test : norme d'un vecteur de dimension 3
	Vector3d v;
	v.setZero();
	v(1) = 2.0;
	EXPECT_EQ(v.rows(), 3);
	EXPECT_EQ(v.cols(), 1);
	EXPECT_DOUBLE_EQ(v(0), 0.0);
	EXPECT_DOUBLE_EQ(v(1), 2.0);
	EXPECT_DOUBLE_EQ(v(2), 0.0);
	EXPECT_DOUBLE_EQ(v.norm(), 2.0);

	// Test : calcul de la norme d'un deuxième vecteur 3D
	Vector3d v2;
	v2(0) = 4.0;
	v2(1) = 2.0;
	v2(2) = 5.0;
	EXPECT_EQ(v2.rows(), 3);
	EXPECT_EQ(v2.cols(), 1);
	EXPECT_DOUBLE_EQ(v2(0), 4.0);
	EXPECT_DOUBLE_EQ(v2(1), 2.0);
	EXPECT_DOUBLE_EQ(v2(2), 5.0);
	EXPECT_DOUBLE_EQ(v2.norm(), 6.7082039324993690892275210061938);

	// Test : produit scalaire 
	EXPECT_DOUBLE_EQ(v.dot(v2), 4.0);

	// Test : matrice identité 4x4
	Matrix4d M;
	M.setIdentity();
	EXPECT_DOUBLE_EQ(M(0, 0), 1.0);
	EXPECT_DOUBLE_EQ(M(0, 1), 0.0);
	EXPECT_DOUBLE_EQ(M(0, 2), 0.0);
	EXPECT_DOUBLE_EQ(M(1, 1), 1.0);
	EXPECT_DOUBLE_EQ(M(1, 0), 0.0);
	EXPECT_DOUBLE_EQ(M(1, 2), 0.0);
	EXPECT_DOUBLE_EQ(M(2, 0), 0.0);
	EXPECT_DOUBLE_EQ(M(2, 1), 0.0);
	EXPECT_DOUBLE_EQ(M(2, 2), 1.0);

	// Test : création d'une matrice de rotation de 45 degrés autour de l'axe des x
	const auto Rx = makeRotation<double>(M_PI / 4.0, 0, 0);
	EXPECT_NEAR(Rx(0, 0), 1, 1e-3); EXPECT_NEAR(Rx(0, 1), 0, 1e-3); EXPECT_NEAR(Rx(0, 2), 0, 1e-3);
	EXPECT_NEAR(Rx(1, 0), 0, 1e-3); EXPECT_NEAR(Rx(1, 1), 0.7071, 1e-3); EXPECT_NEAR(Rx(1, 2), -0.7071, 1e-3);
	EXPECT_NEAR(Rx(2, 0), 0, 1e-3); EXPECT_NEAR(Rx(2, 1), 0.7071, 1e-3); EXPECT_NEAR(Rx(2, 2), 0.7071, 1e-3);

	// Test : création d'une matrice de rotation de 45 degrés autour de l'axe des y
	const auto Ry = makeRotation<double>(0, M_PI / 4.0, 0);
	EXPECT_NEAR(Ry(0, 0), 0.7071, 1e-3); EXPECT_NEAR(Ry(0, 1), 0, 1e-3); EXPECT_NEAR(Ry(0, 2), 0.7071, 1e-3);
	EXPECT_NEAR(Ry(1, 0), 0, 1e-3); EXPECT_NEAR(Ry(1, 1), 1, 1e-3); EXPECT_NEAR(Ry(1, 2), 0, 1e-3);
	EXPECT_NEAR(Ry(2, 0), -0.7071, 1e-3); EXPECT_NEAR(Ry(2, 1), 0, 1e-3); EXPECT_NEAR(Ry(2, 2), 0.7071, 1e-3);

	// Test : création d'une matrice de rotation de 45 degrés autour de l'axe des z
	const auto Rz = makeRotation<double>(0, 0, M_PI / 4.0);
	EXPECT_NEAR(Rz(0, 0), 0.7071, 1e-3); EXPECT_NEAR(Rz(0, 1), -0.7071, 1e-3); EXPECT_NEAR(Rz(0, 2), 0, 1e-3);
	EXPECT_NEAR(Rz(1, 0), 0.7071, 1e-3); EXPECT_NEAR(Rz(1, 1), 0.7071, 1e-3); EXPECT_NEAR(Rz(1, 2), 0, 1e-3);
	EXPECT_NEAR(Rz(2, 0), 0, 1e-3); EXPECT_NEAR(Rz(2, 1), 0, 1e-3); EXPECT_NEAR(Rz(2, 2), 1, 1e-3);

	// Test : création d'une matrice de rotation quelconque.
	const auto Rxyz = makeRotation<double>(M_PI / 3.0, -M_PI / 6.0, M_PI / 4.0);
	EXPECT_NEAR(Rxyz(0, 0), 0.6124, 1e-3); EXPECT_NEAR(Rxyz(0, 1), -0.6597, 1e-3); EXPECT_NEAR(Rxyz(0, 2), 0.4356, 1e-3);
	EXPECT_NEAR(Rxyz(1, 0), 0.6124, 1e-3); EXPECT_NEAR(Rxyz(1, 1), 0.0474, 1e-3); EXPECT_NEAR(Rxyz(1, 2), -0.7891, 1e-3);
	EXPECT_NEAR(Rxyz(2, 0), 0.5, 1e-3); EXPECT_NEAR(Rxyz(2, 1), 0.75, 1e-3); EXPECT_NEAR(Rxyz(2, 2), 0.4330, 1e-3);

	// Test : création d'une transformation homogène via la sous-matrice 3x3 en
	// utilisant la fonction `block`
	M.block(0, 0, 3, 3) = Rxyz;
	M(0, 3) = -0.1;
	M(1, 3) = 1.0;
	M(2, 3) = 2.1;

	// Test : calcule l'inverse de la matrice M et vérifie que M^(-1) * M * v = v
	const Matrix4d Minv = M.inverse();
	const Vector3d v3 = Minv * (M * v2);
	EXPECT_DOUBLE_EQ(v3(0), v2(0));
	EXPECT_DOUBLE_EQ(v3(1), v2(1));
	EXPECT_DOUBLE_EQ(v3(2), v2(2));

	// Test : translation d'un vecteur 3D effectuée avec une matrice 4x4 en coordonnées homogènes
	Matrix4d T;
	T.setIdentity();
	T(0, 3) = 1.2;
	T(1, 3) = 2.5;
	T(2, 3) = -4.0;
	const Vector3d t = T * v3;
	EXPECT_DOUBLE_EQ(t(0), v3(0) + 1.2);
	EXPECT_DOUBLE_EQ(t(1), v3(1) + 2.5);
	EXPECT_DOUBLE_EQ(t(2), v3(2) - 4.0);

	// Test : inverse d'un matrice de rotation
	const Matrix3d Rinv = Rxyz.inverse();
	const Matrix3d RT = Rxyz.transpose<double, 3, 3, ColumnStorage>();
	EXPECT_DOUBLE_EQ(Rinv(0, 0), RT(0, 0));
	EXPECT_DOUBLE_EQ(Rinv(1, 1), RT(1, 1));
	EXPECT_DOUBLE_EQ(Rinv(0, 2), RT(0, 2));


}

/**
 * Test des performance de la multiplication  matrice * vecteur
 * pour de grandes dimensions.
 */

TEST(TestLabo1, PerformanceMatrixVector)
{
	Matrix<double> A(16384, 16384);     // grande matrice avec stockage colonne
	Vector<double> v(16384);            // grand vecteur

	using namespace std::chrono;
	// Test : multiplication avec l'algorithme naif.
	high_resolution_clock::time_point t = high_resolution_clock::now();
	naiveMatrixMult(A, v);
	const duration<double> naive_t = duration_cast<duration<double>>(high_resolution_clock::now() - t);

	// Test : multiplication avec l'implémentation spécifique pour les matrices avec
	// stockage par colonnes.
	t = high_resolution_clock::now();
	A* v;
	const duration<double> optimal_t = duration_cast<duration<double>>(high_resolution_clock::now() - t);

	EXPECT_TRUE(optimal_t < 0.4 * naive_t)
		<< "Naive time: " << duration_cast<std::chrono::milliseconds>(naive_t).count() << " ms, "
		<< "optimized time: " << duration_cast<std::chrono::milliseconds>(optimal_t).count() << " ms";
}

/**
 * Test des performances de l'addition  matrice + matrice
 * pour de grandes dimensions.
 */

TEST(TestLabo1, PerformanceLargeMatrixMatrix)
{
	// deux grandes matrices à stockage par colonnes
	Matrix<double> A(16384, 16384);
	Matrix<double> B(16384, 16384);

	using namespace std::chrono;
	high_resolution_clock::time_point t = high_resolution_clock::now();
	// Test : addition avec l'algorithme naif 
	naiveMatrixAddition(A, B);
	const duration<double> naive_t = duration_cast<duration<double>>(high_resolution_clock::now() - t);

	// Test : addition avec l'implémentation spécifique pour les matrices à
	// stockage par colonnes.
	t = high_resolution_clock::now();
	A + B;
	const duration<double> optimal_t = duration_cast<duration<double>>(high_resolution_clock::now() - t);

	EXPECT_TRUE(optimal_t < 0.4 * naive_t);
}

TEST(TestLabo1, Supplementaires)
{
	// TODO remplacez le code avec vos propres tests ici

	Matrix<double, Dynamic, Dynamic, RowStorage> RowM(10, 20);
	RowM.setZero();
	RowM(0, 0) = 2.1;
	RowM(2, 2) = 9.5;
	RowM(2, 3) = 7.4;
	RowM(3, 2) = 8.5;
	RowM(3, 3) = 6.3;
	RowM(4, 3) = 12.3;

	// Test 1: vérifie le bon fonctionne d'une sous-matrice créer à partir d'une matrice
	SubMatrix<double, Dynamic, Dynamic, RowStorage> SubM(RowM, 2, 2, 3, 2);
	

	EXPECT_EQ(SubM.rows(), 3);
	EXPECT_EQ(SubM.cols(), 2);
	EXPECT_DOUBLE_EQ(SubM(0, 0), RowM(2, 2));
	EXPECT_DOUBLE_EQ(SubM(0, 0), 9.5);
	EXPECT_DOUBLE_EQ(SubM(1, 1), RowM(3, 3));
	EXPECT_DOUBLE_EQ(SubM(1, 1), 6.3);
	EXPECT_DOUBLE_EQ(SubM(2, 1), RowM(4, 3));
	EXPECT_DOUBLE_EQ(SubM(2, 1), 12.3);

	// Test 2: vérifie la copie d'une sous-matrice dans une autre sous-matrice
	SubMatrix<double, Dynamic, Dynamic, RowStorage> SubM2 = SubM;
	EXPECT_EQ(SubM2.rows(), 3);
	EXPECT_EQ(SubM2.cols(), 2);
	EXPECT_DOUBLE_EQ(SubM2(0, 0), SubM(0, 0));
	EXPECT_DOUBLE_EQ(SubM2(0, 0), 9.5);
	EXPECT_DOUBLE_EQ(SubM2(1, 1), SubM(1, 1));
	EXPECT_DOUBLE_EQ(SubM2(1, 1), 6.3);

	
	// Test 3: vérifie la copie d'une matrice dans une sous-matrice
	Matrix<double, Dynamic, Dynamic, RowStorage> RowM2(3, 2);
	RowM2.setZero();
	RowM2(0, 0) = 44.0;
	RowM2(1, 1) = 66.0;
	RowM2(2, 1) = 52.0;

	SubM = RowM2; 

	EXPECT_DOUBLE_EQ(SubM(0, 0), RowM2(0, 0));
	EXPECT_DOUBLE_EQ(SubM(0, 0), 44.0);
	EXPECT_DOUBLE_EQ(SubM(1, 0), 0.0);
	EXPECT_DOUBLE_EQ(SubM(2, 1), RowM2(2, 1));
	EXPECT_DOUBLE_EQ(SubM(2, 1), 52.0);

	// Test 4: vérifie la création d'une correcte matrice transposée à partir d'une sous-matrice
	Matrix<double, Dynamic, Dynamic, RowStorage> SubM2T = SubM2.transpose<double, Dynamic, Dynamic, RowStorage>();
	EXPECT_EQ(SubM2T.rows(), 2);
	EXPECT_EQ(SubM2T.cols(), 3);
	EXPECT_DOUBLE_EQ(SubM2T(0, 0), SubM2(0, 0));
	EXPECT_DOUBLE_EQ(SubM2T(1, 2), SubM2(2, 1));


	// Test 5: vérifie la creation de la transposée à partir d'une matrice avec un stockage en colonne
	Matrix<double, Dynamic, Dynamic, ColumnStorage> ColM(100, 100);
	ColM.setZero();
	ColM(0, 0) = 1.0;
	ColM(10, 5) = 5.0;
	ColM(11, 6) = 6.0;
	ColM(12, 7) = 7.0;
	ColM(12, 5) = 1.0;
	ColM(11, 5) = 4.0;
	ColM(99, 99) = 99.0;
	ColM(10, 33) = 5.0;

	const auto ColMT = ColM.transpose<double, Dynamic, Dynamic, ColumnStorage>();
	EXPECT_DOUBLE_EQ(ColMT(0, 0), ColM(0, 0));
	EXPECT_DOUBLE_EQ(ColMT(33, 10), ColM(10, 33));

	// Test 6: vérifie la creation d'une matrice identité à partir d'une matrice et une sous-matrice
	Matrix< double, Dynamic, Dynamic, ColumnStorage > ColM2I;
	ColM2I = ColM.block(10, 5, 3, 3);
	ColM2I.setIdentity();

	EXPECT_DOUBLE_EQ(ColM2I(0, 0), 1.0);
	EXPECT_DOUBLE_EQ(ColM2I(1, 1), 1.0);
	EXPECT_DOUBLE_EQ(ColM2I(2, 2), 1.0);

	// Test 7: vérifie que la matrice se redimensionne correctement lorsque les matrices ne sont pas de même taille lors d'une copie
	Matrix< double, Dynamic, Dynamic, ColumnStorage > ColM3(10, 10);
	ColM3.setZero();
	ColM3 = ColM.block(10, 5, 3, 3);
	EXPECT_EQ(ColM3.rows(), 3);
	EXPECT_EQ(ColM3.cols(), 3);

	// Test 8: vérifie la création d'un vecteur à partir d'une liste
	Vector<double> v1(4);
	double arr[4] = { 1.5, -2.3, 3.7, -4.1 };
	v1 = arr;
	EXPECT_DOUBLE_EQ(v1(0), arr[0]);
	EXPECT_DOUBLE_EQ(v1(1), arr[1]);
	EXPECT_DOUBLE_EQ(v1(2), arr[2]);
	EXPECT_DOUBLE_EQ(v1(3), arr[3]);

	// Test 9: vérifie que le calcul de la norme de manhattan du vecteur est bien calculé
	EXPECT_DOUBLE_EQ(v1.normManhattan(), 11.6);

	// Test 10: vérifie que le calcul de la norme infinie et de sa norme-p du vecteur est bien calculé
	EXPECT_DOUBLE_EQ(v1.normInfinie(), 4.1);
	EXPECT_DOUBLE_EQ(v1.normP(3), 5.1313967320559923);

	// Test 11: vérifie la copie d'un vecteur à partir d'une sous matrice
	Matrix< double, Dynamic, Dynamic, ColumnStorage > m1(2, 4);
	m1.setZero();
	m1(1, 1) = 12.5;
	m1(1, 2) = 78.2;

	Vector<double> v2(4);
	v2 = m1.block(0, 1, 4, 1);
	EXPECT_DOUBLE_EQ(v2(0), 0.0);
	EXPECT_DOUBLE_EQ(v2(1), 12.5);

	// Test 12: vérifie la soustraction entre vecteur
	Vector<double> v3(4);
	v3 = v1 - v2;
	EXPECT_DOUBLE_EQ(v3(0), 1.5);
	EXPECT_DOUBLE_EQ(v3(1), -14.8);

	// Test 13: vérifie la soustraction entre matrice générique
	Matrix<double> m2(6, 6);
	m2(4, 3) = 8.5;
	m2(5, 2) = 4.5;
	const auto m22 = m2 - (2.0 * m2);
	EXPECT_DOUBLE_EQ(m22(4, 3), -8.5);
	EXPECT_DOUBLE_EQ(m22(5, 2), -4.5);

	// Test 14: vérifie la soustraction entre matrice avec un stockage en ligne
	Matrix<double, Dynamic, Dynamic, RowStorage> m3(5, 5);
	m3(1, 2) = 5.5;
	m3(2, 3) = 6.0;
	const auto m4 = m3 - (2.0 * m3);
	EXPECT_DOUBLE_EQ(m4(1, 2), -5.5);
	EXPECT_DOUBLE_EQ(m4(2, 3), -6.0);

	// Test 15: vérifie la soustraction entre matrice avec un stockage en colonne
	Matrix<double, Dynamic, Dynamic, ColumnStorage> m5 = m3.transpose();
	const auto m6 = m5 - (2.0 * m5);
	EXPECT_DOUBLE_EQ(m6(2, 1), -5.5);
	EXPECT_DOUBLE_EQ(m6(3, 2), -6.0);

	// Test 17: vérifie la création d'une matrice de mise en échelle
	Matrix<double, 3, 3, ColumnStorage> mScaling(3, 3);
	mScaling = scaling(3.0);
	EXPECT_DOUBLE_EQ(mScaling(0, 0), 3.0);
	EXPECT_DOUBLE_EQ(mScaling(1, 1), 3.0);
	EXPECT_DOUBLE_EQ(mScaling(2, 2), 3.0);


	// Test 18: vérifie la mise en échelle d'une matrice
	Matrix<double, 3, 3, ColumnStorage> m7(3, 3);
	m7(2, 1) = 4;
	m7(1, 1) = 6;
	const auto m7Scale = m7 * mScaling;
	EXPECT_DOUBLE_EQ(m7Scale(2, 1), 12);
	EXPECT_DOUBLE_EQ(m7Scale(1, 1), 18);

	// Test 19 : vérifie le produit vectoriel de deux vecteur 3d
	Vector<double, 3> vi;
	Vector<double, 3> vj;
	double iArr[3] = { 1.0, 0.0, 0.0 };
	double jArr[3] = { 0.0, 1.0, 0.0 };
	vi = iArr;
	vj = jArr;

	const auto vk = cross(vi, vj);
	EXPECT_DOUBLE_EQ(vk(0), 0.0);
	EXPECT_DOUBLE_EQ(vk(1), 0.0);
	EXPECT_DOUBLE_EQ(vk(2), 1.0);

	// Test 20: vérifie si le résultat du produit vectoriel est bien perpendiculaire aux deux autres vecteurs
	EXPECT_DOUBLE_EQ(vi.dot(vj), 0.0);
	EXPECT_DOUBLE_EQ(vi.dot(vk), 0.0);
	EXPECT_DOUBLE_EQ(vj.dot(vk), 0.0);

}

int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	const int ret = RUN_ALL_TESTS();

	return ret;
}
