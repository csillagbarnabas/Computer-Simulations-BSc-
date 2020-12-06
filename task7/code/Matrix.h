#include <vector>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <initializer_list>
#include <cmath>
#include <ostream>
#include <sstream>

struct Idx1{};
struct Idx2{};

namespace detail
{
	template<typename V1, typename V2, typename F>
	void transform_matrix1(V1 const& v1, V2& v2, F f)
	{
		std::transform(v1.cbegin(), v1.cend(), v2.begin(), f);
	}

	template<typename V1, typename V2, typename V3, typename F>
	void transform_matrix2(V1 const& v1, V2 const& v2, V3& v3, F f)
	{
		std::transform(v1.cbegin(), v1.cend(), v2.cbegin(), v3.begin(), f);
	}
}

auto add = [](auto const& x, auto const& y){ return x + y; };
auto sub = [](auto const& x, auto const& y){ return x - y; };

template<typename T>
class Matrix{
	int N;
    std::vector<T> data;
	public:
    	T&  operator()(int i, int j)
    	{ return data[N*i+j]; }
    	T const& operator()(int i, int j) const
    	{ return data[N*i+j]; }
		T&  operator[](int i)
    	{ return data[i]; }
    	T const& operator[](int i) const
    	{ return data[i]; }
	Matrix(): N(0), data(0) {};
	Matrix( Matrix const& ) = default;
	Matrix(Matrix&& m) : N{m.N}, data{std::move(m.data)} {m.N = 0;  };  
	Matrix<T>& operator=(Matrix const&) = default;
	Matrix<T>& operator=(Matrix && m){
		if(&data == &m.data){
			return *this;
		}
		else{
			N=m.N;
			data=std::move(m.data);
			m.N = 0; 
			return *this;
		}
	}
	template<typename F>
	Matrix(Idx1, F f,int M){
		data.resize(M*M);
		N = M;
		for(int i=0; i<M*M; ++i){
			data[i] = f(i);
		}
	}
	template<typename F>
	Matrix(Idx2, F f, int M){
		data.resize(M*M);
		N = M;
		for(int i=0; i<M; ++i){
			for( int j=0; j < M; ++j){
				data[i*M+j] = f(i,j);
			}
		}
	}
	Matrix( int n, std::vector<T> const& x ) : N(n), data(x){};
    Matrix<T>& operator+= (Matrix<T> const& cpy)
	{
		detail::transform_matrix2(data, cpy.data, data, add);
		return *this;
	}
	Matrix<T>& operator-= (Matrix<T> const& cpy)
	{
		detail::transform_matrix2(data, cpy.data, data, sub);
		return *this;
	}
	Matrix<T>& operator*= (T const& scl)
	{
		detail::transform_matrix1(data, data, [=](T const& x){ return x * scl;} );
		return *this;
	}
	Matrix<T>& operator/= (T const& scl)
	{
		detail::transform_matrix1(data, data, [=](T const& x){ return x / scl;} );
		return *this;
	}
	int datasize()const{
		return static_cast<int>(data.size());
		}
	int Nsize()const{
		return N;
		}
	auto begin()
	{
		return data.begin();
	}

	auto cbegin() const
	{
		return data.cbegin();
	}

	auto end()
	{
		return data.end();
	}

	auto cend() const
	{
		return data.cend();
	}
    };
template<typename T>
Matrix<T> operator+(Matrix<T> const& A, Matrix<T> const& B)
{
	return Matrix<T>(Idx1{}, [&](auto i){ return A[i] + B[i]; },A.Nsize());
}
template<typename T>
Matrix<T> && operator+( Matrix<T> && m1, Matrix<T> const& m2 )
{
	detail::transform_matrix2(m1, m2, m1, add);
	return std::move(m1);
}
template<typename T>
Matrix<T> && operator+( Matrix<T> const& m1, Matrix<T> && m2 )
{
	detail::transform_matrix2(m1, m2, m2, add);
	return std::move(m2);
}
template<typename T>
Matrix<T> && operator+( Matrix<T> && m1, Matrix<T> && m2 )
{
	detail::transform_matrix2(m1, m2, m1, add);
	return std::move(m1);
}
//--
template<typename T>
Matrix<T> operator-(Matrix<T> const& A, Matrix<T> const& B)
{
	return Matrix<T>(Idx1{}, [&](auto i){ return A[i] - B[i]; },A.Nsize());
}
template<typename T>
Matrix<T> && operator-( Matrix<T> && m1, Matrix<T> const& m2 )
{
	detail::transform_matrix2(m1, m2, m1, sub);
	return std::move(m1);
}
template<typename T>
Matrix<T> && operator-( Matrix<T> const& m1, Matrix<T> && m2 )
{
	detail::transform_matrix2(m1, m2, m2, sub);
	return std::move(m2);
}
template<typename T>
Matrix<T> && operator-( Matrix<T> && m1, Matrix<T> && m2 )
{
	detail::transform_matrix2(m1, m2, m1, sub);
	return std::move(m1);
}
template<typename T>
Matrix<T>&& operator*(Matrix<T> && m, T const& scl)
{
	detail::transform_matrix1(m, m, [=](T const& x){ return x * scl; });
	return std::move(m);
}
template<typename T>
Matrix<T> operator*(Matrix<T> const& A, T const& s)
{
	return Matrix<T>(Idx1{}, [&](auto i){ return A[i] * s; },A.Nsize());
}
template<typename T>
Matrix<T> operator*(T const& s,Matrix<T> const& A)
{
	return Matrix<T>(Idx1{}, [&](auto i){ return s * A[i]; },A.Nsize());
}
template<typename T>
Matrix<T>&& operator*(T const& scl,Matrix<T> && m)
{
	detail::transform_matrix1(m, m, [=](T const& x){ return scl * x; });
	return std::move(m);
}
template<typename T>
Matrix<T> operator/(Matrix<T> const& A, T const& s)
{
	return Matrix<T>(Idx1{}, [&](auto i){ return A[i] / s; },A.Nsize());
}
template<typename T>
Matrix<T>&& operator/(Matrix<T> && m, T const& scl)
{
	detail::transform_matrix1(m, m, [=](T const& x){ return x / scl; });
	return std::move(m);
}
template<typename T>
Matrix<T> operator*(Matrix<T> const& A, Matrix<T> const& B)
{
	auto N = A.Nsize();
	return Matrix<T>(Idx2{}, [&](int i, int j)
	{
	T sum = 0.0;
	for(int k = 0; k<N; ++k){ sum += A(i, k) * B(k, j); }
	return sum;
	}, N);
}
template<typename T>
Matrix<T>&& operator*(Matrix<T> && A, Matrix<T> const& B)
{
	std::vector<T> v(A.Nsize());
	for(int i=0; i<B.Nsize(); i++){
		for(int j=0; j<B.Nsize(); j++){
			T h = 0.0;
			for(int k=0; k<B.Nsize(); k++){
				h += A(i,k) * B(k,j);
			}
			v[j] = h;
		}
		for(int j=0; j<A.Nsize(); j++){
			A(i,j) = v[j];
		}
	}
	return std::move(A);
}
template<typename T>
Matrix<T>&& operator*(Matrix<T> const& A, Matrix<T> && B)
{
	std::vector<T> v(A.Nsize());
	for(int j=0; j<A.Nsize(); j++){
		for(int i=0; i<A.Nsize(); i++){
			T h = 0.0;
			for(int k=0; k<A.Nsize(); k++){
				h += A(j,k) * B(k,i);
			}
			v[i] = h;
		}
		for(int i=0; i<B.Nsize(); i++){
			B(j,i) = v[i];
		}
	}
	return std::move(B);
}
template<typename T>
Matrix<T>&& operator*(Matrix<T> && A, Matrix<T> && B)
{
	std::vector<T> v(A.Nsize());
	for(int i=0; i<B.Nsize(); i++){
		for(int j=0; j<A.Nsize(); j++){
			T h = 0.0;
			for(int k=0; k<A.Nsize(); k++){
				h += A(i,k) * B(k,j);
			}
			v[j] = h;
		}
		for(int j=0; j<A.Nsize(); j++){
			A(i,j) = v[j];
		}
	}
	return std::move(A);
}
std::istream& operator>>( std::istream& s, Matrix<double>& A )
{
	std::string tmp;
	for(int j=0; j<A.Nsize(); j++){
		std::getline(s, tmp);
		if(tmp.size() > 0)
		{
		std::stringstream ss(tmp);
		for(int k=0; k<A.Nsize(); k++){
			if(k < A.Nsize()-1){
				std::getline(ss, tmp, ',');
				A(j,k) = std::stod(tmp);
			}
			else{
				std::getline(ss, tmp, '\n');
				A(j,k) = std::stod(tmp);
				}
			}
		}
	}
	return s;
}
template<typename T>
std::ostream& operator<< (std::ostream& o, Matrix<T> const& m)
{
	for(int j = 0; j < m.Nsize(); ++j){
		for(int i=0; i<m.Nsize(); ++i)
		{
			if(m(j,i) >= 0){
				if(i != m.Nsize()-1){
					o << " " << m(j,i) << ",";
				}
				else{
					o << " " << m(j,i) << std::endl;
				}
			}
			if(m(j,i) < 0){
				if(i != m.Nsize()-1){
					o << m(j,i) << ",";
				}
				else{
					o << m(j,i) << std::endl;
				}
			}
		}
	}
	return o;
}
