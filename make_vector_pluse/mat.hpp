#ifndef KMUVCL_GRAPHICS_MAT_HPP
#define KMUVCL_GRAPHICS_MAT_HPP

#include <iostream>
#include <cstring>
#include <cstdarg>

namespace kmuvcl {
  namespace math {

    template <unsigned int M, unsigned int N, typename T>
    class mat
    {
    public:
      mat()
      {
        set_to_zero();
      }

      mat(const T elem)
      {
        std::fill(val, val + M*N, elem);
      }

      T& operator()(unsigned int r, unsigned int c)
      {
        return val[r+c*M];
      }

      const T& operator()(unsigned int r, unsigned int c) const
      {
        return val[r+c*M];
      }

      // type casting operators
      operator const T* () const
      {
        return  val;
      }

      operator T* ()
      {
        return  val;
      }

      void set_to_zero()
      {
        int MaxIndex = M*N;
        for(int i=0; i< MaxIndex; ++i)
        {
          val[i] = 0;
        }
      }

      void get_ith_column(unsigned int i, vec<M, T>& col) const
      {
        int k = 0;
        int startIndex = i*M;
        int endIndex = i*M+M;
    
        for(int j=startIndex; j < endIndex; ++j)
        {
          col[k] = val[j]; 
          ++k;
        }
      }

      void set_ith_column(unsigned int i, const vec<M, T>& col)
      {
        int k = 0;
        int startIndex = i*M;
        int endIndex = i*M+M;
    
        for(int j=startIndex; j < endIndex; ++j)
        {
          val[j] = col[k]; 
          ++k;
        }
      }

      void get_ith_row(unsigned int i, vec<N, T>& row) const
      {
        int k = 0;
        int startIndex = i;
        int endIndex = i+N*N;
        
        for(int j=startIndex; j < endIndex; j+=N)
        {
          row[k] = val[j]; 
          ++k;
        }
      }

      void set_ith_row(unsigned int i, const vec<N, T>& row)
      {
        int k = 0;
        int startIndex = i;
        int endIndex = i+N*N;
        
        for(int j=startIndex; j < endIndex; j+=N)
        {
          val[j] = row[k]; 
          ++k;
        }
      }

      mat<N, M, T> transpose() const
      {
        mat<N, M, T>  trans;
        for(int i=0; i<N;++i)
        {
          for(int j=0; j< M; ++j)
          {
            trans(i,j) = this->val[i*M+j];
          }
        }        
        return  trans;
      }

    protected:
      T val[M*N];   // column major
    };

    typedef mat<3, 3, float>    mat3x3f;
    typedef mat<3, 3, double>   mat3x3d;

    typedef mat<4, 4, float>    mat4x4f;
    typedef mat<4, 4, double>   mat4x4d;

  } // math
} // kmuvcl

#include "operator.hpp"

#endif // #ifndef KMUVCL_GRAPHICS_MAT_HPP
