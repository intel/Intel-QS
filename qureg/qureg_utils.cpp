//------------------------------------------------------------------------------
// Copyright 2017 Intel Corporation
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//------------------------------------------------------------------------------

#include "qureg.hpp"


/// \addtogroup qureg
/// @{

/// @file qureg_utils.cpp
///  @brief Define the @c QubitRegister methods used as basic operations.

/////////////////////////////////////////////////////////////////////////////////////////
/// @brief ???
template <class Type>
bool QubitRegister<Type>::operator==(const QubitRegister &rhs)
{
  assert(rhs.GlobalSize() == GlobalSize());
  for (std::size_t i = 0; i < rhs.LocalSize(); i++)
  {
      if (state[i] != rhs.state[i])
      {
          printf("(%lf,%lf) != (%lf,%lf)\n", state[i].real(), state[i].imag(),
                                             rhs.state[i].real(), rhs.state[i].imag());
          return false;
      }
  }
  return true;
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief ???
template <class Type>
QubitRegister<Type>::BaseType QubitRegister<Type>::maxabsdiff(QubitRegister &x, Type sfactor)
{
  assert(LocalSize() == x.LocalSize());
  BaseType lcl_maxabsdiff = -1.0;

  std::size_t lcl = LocalSize();
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for simd reduction(max : lcl_maxabsdiff)
#else
   TODO(Remember to find 'omp parallel for simd reduction' equivalent for gcc)
#endif
  for (std::size_t i = 0; i < lcl; i++) {
    lcl_maxabsdiff = std::max(lcl_maxabsdiff, std::abs(state[i] - sfactor*x.state[i]));
  }

  BaseType glb_maxabsdiff ;
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = openqu::mpi::Environment::comm();
  // MPI_Allreduce(&lcl_maxabsdiff, &glb_maxabsdiff, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce_x(&lcl_maxabsdiff, &glb_maxabsdiff,  MPI_MAX, comm);
#else
  glb_maxabsdiff = lcl_maxabsdiff;
#endif

  return glb_maxabsdiff;
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief ???
template <class Type>
QubitRegister<Type>::BaseType QubitRegister<Type>::maxl2normdiff(QubitRegister &x)
{
  assert(LocalSize() == x.LocalSize());
  BaseType lcl_diff = 0.0;
  // #pragma omp parallel for simd reduction(+:lcl_diff)
  std::size_t lcl = LocalSize();
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for reduction(+ : lcl_diff)
#else
   TODO(Remember to find 'omp parallel for simd reduction' equivalent for gcc)
#endif
  for (std::size_t i = 0; i < lcl; i++)
  {
      Type r = state[i] - x.state[i];
      lcl_diff += std::norm(r);
  }

  BaseType glb_diff;
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = openqu::mpi::Environment::comm();
  // MPI_Allreduce(&lcl_diff, &glb_diff, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce_x(&lcl_diff, &glb_diff,  MPI_MAX, comm);
#else
   glb_diff = lcl_diff;
#endif

  return glb_diff;
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Normalize the quantum state (L2 norm).
template <class Type>
void QubitRegister<Type>::Normalize() 
{
  BaseType global_norm = ComputeNorm();
  std::size_t lcl = LocalSize();
#pragma omp parallel for 
  for(std::size_t i = 0; i < lcl; i++)
  {
     state[i] = state[i] / global_norm;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Compute the norm of the state (L2 norm).
template <class Type>
QubitRegister<Type>::BaseType QubitRegister<Type>::ComputeNorm()
{
  BaseType local_normsq = 0;
  std::size_t lcl = LocalSize();
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for reduction(+ : local_normsq)
#else
   TODO(Remember to find 'omp parallel for simd reduction' equivalent for gcc)
#endif
  for(std::size_t i = 0; i < lcl; i++)
  {
     local_normsq += std::norm(state[i]);
  }

  BaseType global_normsq;
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = openqu::mpi::Environment::comm();
  // MPI_Allreduce(&local_normsq, &global_normsq, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce_x(&local_normsq, &global_normsq,  MPI_SUM, comm);
#else
  global_normsq = local_normsq;
#endif

  return std::sqrt(global_normsq);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Compute the overlap <psi|this state>
/// @param psi is the second state
///
/// The overlap between this state and another state |psi\>
/// is define by:\n
///     \<psi|this state\>
template <class Type>
Type QubitRegister<Type>::ComputeOverlap( QubitRegister<Type> &psi)
{
  Type local_over = Type(0.,0.);
  BaseType local_over_re = 0.;
  BaseType local_over_im = 0.;
  std::size_t lcl = LocalSize();
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for private(local_over) reduction(+ : local_over_re,local_over_im)
#else
   TODO(Remember to find 'omp parallel for simd reduction' equivalent for gcc)
#endif
  for(std::size_t i = 0; i < lcl; i++)
  {
     local_over = std::conj(psi[i]) * state[i] ; 
     local_over_re +=  std::real( local_over );
     local_over_im +=  std::imag( local_over );
  }
  
  BaseType global_over_re(0.) , global_over_im(0.) ;
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = openqu::mpi::Environment::comm();
  // MPI_Allreduce(&local_over_re, &global_over_re, 1, MPI_DOUBLE, MPI_SUM, comm);
  // MPI_Allreduce(&local_over_im, &global_over_im, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce_x(&local_over_re, &global_over_re,  MPI_SUM, comm);
  MPI_Allreduce_x(&local_over_im, &global_over_im,  MPI_SUM, comm);
#else
  global_over_re = local_over_re;
  global_over_im = local_over_im;
#endif

  return Type(global_over_re,global_over_im);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief ???
template <class Type>
double QubitRegister<Type>::Entropy()
{

  std::size_t lcl = LocalSize();
  double local_Hp = 0;

  if(timer) timer->Start("ENT", 0);

  double ttot = 0., ttmp1 = sec();
#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for reduction(+ : local_Hp)
#else
   TODO(Remember to find 'omp parallel for simd reduction' equivalent for gcc)
#endif
  for (std::size_t i = 0; i < lcl; i++)
  {
      double pj = std::norm(state[i]) ;
      if (pj != double(0.))
      {
          local_Hp -= pj * std::log(pj);
      }
  }

  double global_Hp;
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = openqu::mpi::Environment::comm();
  // MPI_Allreduce(&local_Hp, &global_Hp, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce_x(&local_Hp, &global_Hp,  MPI_SUM, comm);
#else
  global_Hp = local_Hp;
#endif
  
  ttot = sec() - ttmp1;
 
  if (timer)
  {
      double datab = D(sizeof(state[0])) * D(lcl) / ttot;
      timer->record_sn(ttot, datab / ttot);
      timer->Stop();
  }

  return global_Hp / (double)log(double(2.0));
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief ???
template <class Type>
std::vector<double> QubitRegister<Type>::GoogleStats()
{
  std::vector <double> stats;

  std::size_t lcl = LocalSize();
  double two2n = D(GlobalSize());
  
  double entropy = 0, avgselfinfo=0,
         m2 = 0, m3 = 0, m4 = 0, m5 = 0, m6 = 0, 
         m7 = 0, m8 = 0, m9 = 0, m10 = 0; 

  if(timer) timer->Start("ENT", 0);

  double ttot = 0., ttmp1 = sec();

#if defined(__ICC) || defined(__INTEL_COMPILER)
#pragma omp parallel for reduction(+ : entropy, avgselfinfo, m2, m3, m4, m5, m6, m7, m8, m9, m10)
#else
   TODO(Remember to find 'omp parallel for simd reduction' equivalent for gcc)
#endif
  #pragma simd reduction(+ : entropy, avgselfinfo, m2, m3, m4, m5, m6, m7, m8, m9, m10)
  // #pragma novector
  for (std::size_t i = 0; i < lcl; i++)
  {
    double pj = std::norm(state[i]) ;
    if (pj != double(0.))
    {
        double nl = log(pj);
        // double nl = pj*pj;
        entropy -= pj * nl;
        avgselfinfo -= nl;
    }
    double pj2  = pj *  pj,
           pj3  = pj2 * pj,
           pj4  = pj2 * pj2,
           pj5  = pj3 * pj2,
           pj6  = pj3 * pj3,
           pj7  = pj4 * pj3,
           pj8  = pj4 * pj4,
           pj9  = pj5 * pj4,
           pj10 = pj5 * pj5;
    m2  += pj2;
    m3  += pj3;
    m4  += pj4;
    m5  += pj5;
    m6  += pj6;
    m7  += pj7;
    m8  += pj8;
    m9  += pj9;
    m10 += pj10;
  }

  double global_entropy;
  double global_avgselfinfo;
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = openqu::mpi::Environment::comm();
  MPI_Allreduce_x(&entropy, &global_entropy,  MPI_SUM, comm);
  MPI_Allreduce_x(&avgselfinfo, &global_avgselfinfo,  MPI_SUM, comm);
#else
  global_entropy = entropy;
  global_avgselfinfo = avgselfinfo;
#endif
  global_entropy /= (double)std::log(double(2.0));
  stats.push_back(global_entropy);
  global_avgselfinfo /= (double)log(double(2.0));
  global_avgselfinfo /= two2n;
  stats.push_back(global_avgselfinfo);

  // compute moments
  std::vector <double> m = {m2, m3, m4, m5, m6, m7, m8, m9, m10},
                       factor(m.size()), 
                       global_m(m.size());
  double factorial = 1.0;
  for(auto i = 0; i < m.size(); i++)
  {
      auto k = i + 2;
      factorial *= D(k);
      factor[i] = pow(two2n, D(k - 1)) / factorial;

      m[i] *= factor[i];
#ifdef INTELQS_HAS_MPI
      MPI_Allreduce_x(&(m[i]), &(global_m[i]),  MPI_SUM, comm);
#else
      global_m[i] = m[i];
#endif
      stats.push_back(global_m[i]);
  }

  ttot = sec() - ttmp1;

  if (timer)
  {
      double datab = D(sizeof(state[0])) * D(lcl) / ttot;
      timer->record_sn(ttot, datab / ttot);
      timer->Stop();
  }

  return stats;
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Prepare information about the state amplitudes.
/// @param pcum is the cumulative probability for the (possibly partial) state
///
/// Info provided in the output string:
/// index of the computational state and corresponding amplitude.
/// Partial sum of |amplitude|^2 is computed.
//------------------------------------------------------------------------------
template <class Type, class BaseType>
std::string PrintVector(Type *state, std::size_t size, std::size_t num_qubits,
                        BaseType &cumulative_probability, Permutation *permutation)
{
  std::string str;
  int rank = 0;
  for (std::size_t i = 0; i < size; i++) {
    // std::string bin = dec2bin(rank * size + i, num_qubits, false);
    std::string bin = permutation->lin2perm(rank * size + i);
    char s[4096];
    sprintf(s, "\t%-13.8lf + i * %-13.8lf   %% |%s> p=%lf\n",
            std::real(state[i]), std::imag(state[i]),
            bin.c_str(), std::norm(state[i]) );
    str = str + s;
    cumulative_probability += std::norm(state[i]);
  }
  return std::string(str);
}


/// @brief Print on screen some information about the state.
/// @param x is the message to describe the state to be printed
///
/// The pieces of information that are printed are:
/// - permutation map
/// - all amplitudes of the computational basis
/// - cumulative_probability corresponds to the state norm
//------------------------------------------------------------------------------
template <class Type>
void QubitRegister<Type>::Print(std::string x, std::vector<std::size_t> qubits)
{
  TODO(Second argument of Print() is not used!)
  BaseType cumulative_probability = 0;

  unsigned myrank=0, nprocs=1;
#ifdef INTELQS_HAS_MPI
  myrank = openqu::mpi::Environment::rank();
  nprocs = openqu::mpi::Environment::size();
  MPI_Comm comm = openqu::mpi::Environment::comm();
  openqu::mpi::barrier();
#endif

  if (myrank == 0)
  {
      // print permutation
      assert(permutation);
      printf("permutation: %s\n", permutation->GetMapStr().c_str());
      std::string s = PrintVector<Type, BaseType>(state, LocalSize(), num_qubits,
                                                  cumulative_probability, permutation);
      printf("%s=[\n", x.c_str());
      printf("%s", s.c_str());
#ifdef INTELQS_HAS_MPI
      for (std::size_t i = 1; i < nprocs; i++)
      {
          std::size_t len;
#ifdef BIGMPI
          MPIX_Recv_x(&len, 1, MPI_LONG, i, 1000 + i, comm, MPI_STATUS_IGNORE);
#else
          MPI_Recv(&len, 1, MPI_LONG, i, 1000 + i, comm, MPI_STATUS_IGNORE);
#endif //BIGMPI
          s.resize(len);
#ifdef BIGMPI
          MPIX_Recv_x((void *)(s.c_str()), len, MPI_CHAR, i, i, comm, MPI_STATUS_IGNORE);
#else
          MPI_Recv((void *)(s.c_str()), len, MPI_CHAR, i, i, comm, MPI_STATUS_IGNORE);
#endif //BIGMPI
          printf("%s", s.c_str());
      }
#endif
  }
  else
  {
#ifdef INTELQS_HAS_MPI
      std::string s = PrintVector(state, LocalSize(), num_qubits, cumulative_probability, permutation);
      std::size_t len = s.length() + 1;
#ifdef BIGMPI
      MPIX_Send_x(&len, 1, MPI_LONG, 0, 1000 + myrank, comm);
      MPIX_Send_x(const_cast<char *>(s.c_str()), len, MPI_CHAR, 0, myrank, comm);
#else
      MPI_Send(&len, 1, MPI_LONG, 0, 1000 + myrank, comm);
      MPI_Send(const_cast<char *>(s.c_str()), len, MPI_CHAR, 0, myrank, comm);
#endif //BIGMPI
#endif
  }

  BaseType glb_cumulative_probability;
#ifdef INTELQS_HAS_MPI
  MPI_Reduce(&cumulative_probability, &glb_cumulative_probability, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
#else
  glb_cumulative_probability = cumulative_probability;
#endif
  if (myrank == 0)
  {
      printf("]; %% cumulative probability = %lf\n", glb_cumulative_probability);
  }

#ifdef INTELQS_HAS_MPI
  openqu::mpi::barrier();
#endif
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Only used for the MPI implementation.
//------------------------------------------------------------------------------
template <class Type>
void QubitRegister<Type>::dumpbin(std::string fn)
{
#ifdef INTELQS_HAS_MPI
  MPI_Comm comm = openqu::mpi::Environment::comm();
  unsigned myrank = openqu::mpi::Environment::rank();
  MPI_Status status;
  MPI_File fh;
  MPI_File_open(comm, (char *)fn.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &fh);
  std::size_t size = LocalSize();
  assert(size < INT_MAX);
  MPI_Offset offset = size * UL(myrank * sizeof(Type));

  double t0 = sec();
  openqu::mpi::barrier();
  for (int i = 0; i < num_block; i++){
    OutputLclState(i);
    MPI_File_write_at(fh, offset + i * block_size * sizeof(Type), (void *)(&(state[0])), block_size, MPI_DOUBLE_COMPLEX, &status);
  }
  openqu::mpi::barrier();
  double t1 = sec();
  MPI_File_close(&fh);
  if (myrank == 0)
  {
      double bw = D(UL(sizeof(state[0])) * size) / (t1 - t0) / 1e6;
      printf("Dumping state to %s took %lf sec (%lf MB/s)\n", fn.c_str(), t1 - t0, bw);
  }
#else
  assert(0);
#endif
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Enable the collection of statistics (time used in computation and bandwidth).
//------------------------------------------------------------------------------
template <class Type>
void QubitRegister<Type>::EnableStatistics()
{
  unsigned myrank=0, nprocs=1;
#ifdef INTELQS_HAS_MPI
  myrank = openqu::mpi::Environment::rank();
  nprocs = openqu::mpi::Environment::size();
#endif

  assert(timer == NULL);
  timer = new Timer(num_qubits, myrank, nprocs);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Print the statistics (time used in computation and bandwidth) to screen.
//------------------------------------------------------------------------------
template <class Type>
void QubitRegister<Type>::GetStatistics()
{
  assert(timer);
  timer->Breakdown();
  // delete timer;
  // timer = NULL;
}


/////////////////////////////////////////////////////////////////////////////////////////
/// @brief Reset the statistics to allow a new start.
//------------------------------------------------------------------------------
template <class Type>
void QubitRegister<Type>::ResetStatistics()
{
// FIXME it does not delete the timer object!
  assert(timer);
  delete timer;
  timer = NULL;
}


template class QubitRegister<ComplexSP>;
template class QubitRegister<ComplexDP>;

/// @}

int sz_mode = PW_REL;
size_t sz_thres = 10;
double eb_sz = 1E-1;

template <class Type>
void QubitRegister<Type>::OutputStatesToFile(int count)
{
  string filename = "state-" + std::to_string(count) + ".bin";
  dumpbin(filename);
}


#define LEVEL 6
#define THRESHOLD_RATIO 64
int mode[LEVEL] = {0, PW_REL, PW_REL, PW_REL, PW_REL, PW_REL};
double eb[LEVEL] = {0, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1};


template <class Type>
unsigned char* QubitRegister<Type>::BlkCompress(double* data_blk, size_t* comp_len)
{
  unsigned char* comp_blk;
  unsigned char* comp_byte;

  if (compress_level == 0){
    *comp_len = sz_lossless_compress(ZSTD_COMPRESSOR, 1, (unsigned char*)data_blk, (unsigned long) block_size * 16, &comp_byte);
  } else {
    SZ_Init(NULL);
    lossy = 1;
    sz_mode = mode[compress_level];
    eb_sz = eb[compress_level];
    sz_thres = 0;
    
    comp_byte = SZ_compress_qcsim(SZ_DOUBLE, data_blk, comp_len, sz_mode, eb_sz, eb_sz, eb_sz, block_size, 0, 16384, 2);

    SZ_Finalize();
  }
  *comp_len = 1 + (*comp_len);
  comp_blk = (unsigned char*) malloc(*comp_len);
  memcpy(comp_blk, comp_byte, *comp_len);
  free(comp_byte);
  if (compress_level == 0){
    comp_blk[(*comp_len) - 1] = 'Z';
  } else {
    comp_blk[(*comp_len) - 1] = 'S';
  }
  return comp_blk;
}


template <class Type>
void QubitRegister<Type>::BlkDecompress(unsigned char* comp_blk, double* data_blk, size_t comp_len)
{
  if (comp_blk[comp_len - 1] == 'Z'){
    unsigned char * data_byte;
    sz_lossless_decompress(ZSTD_COMPRESSOR, comp_blk, comp_len - 1, &data_byte, (unsigned long)block_size * 16);
    memcpy(data_blk, data_byte, block_size * 16);
    free(data_byte);
  } else {
    SZ_decompress_qcsim(SZ_DOUBLE, comp_blk, comp_len - 1, data_blk, block_size, 2);
  }
  
}

#define SZ_DEBUG 0
#if SZ_DEBUG
int g_count = 0;
#endif

template <class Type>
void QubitRegister<Type>::CompressLclState(int idx)
{
  double t1, t2;
  unsigned rank = openqu::mpi::Environment::rank();
  unsigned nprocs = openqu::mpi::Environment::size();
  MPI_Comm comm = openqu::mpi::Environment::comm();

#if SZ_DEBUG
  char filename[128];
  if (!rank && !idx){
    sprintf(filename, "data-%u-%d-%d", rank, idx, g_count);
    FILE *pFile = fopen(filename, "wb");
    if (pFile == NULL){
      cout << "Failed to open the file: " << filename << endl << flush;
    } else {
      fwrite(state, 1, block_size * 16, pFile);
      fclose(pFile);
    }
  }
#endif

  t1 = MPI_Wtime();
  if (list_compressed_blk[idx] != NULL) free(list_compressed_blk[idx]);
  list_compressed_blk[idx] = BlkCompress((double*) state, &list_len[idx]);

#if SZ_DEBUG
  if (!rank && !idx){
    sprintf(filename, "szcomp-%u-%d-%d", rank, idx, g_count++);
    FILE *pFile = fopen(filename, "wb");
    if (pFile == NULL){
      cout << "Failed to open the file: " << filename << endl << flush;
    } else {
      fwrite(list_compressed_blk[idx], 1, list_len[idx] - 1, pFile);
      fclose(pFile);
    }
  }
#endif

  t2 = MPI_Wtime();
}

template <class Type>
void QubitRegister<Type>::CompressTmpState()
{
  double t1, t2;
  unsigned rank = openqu::mpi::Environment::rank();
  unsigned nprocs = openqu::mpi::Environment::size();
  MPI_Comm comm = openqu::mpi::Environment::comm();

  t1 = MPI_Wtime();
  if (tmp_compressed_blk != NULL) free(tmp_compressed_blk);
  tmp_compressed_blk = BlkCompress((double*) tmp_state, &tmp_len);
  
  t2 = MPI_Wtime();

}

template <class Type>
void QubitRegister<Type>::DecompressLclState(int idx)
{
  double t1, t2;
  int i;
  unsigned rank = openqu::mpi::Environment::rank();
  unsigned nprocs = openqu::mpi::Environment::size();
  MPI_Comm comm = openqu::mpi::Environment::comm();

  t1 = MPI_Wtime();

  BlkDecompress(list_compressed_blk[idx], (double*) state, list_len[idx]);
  t2 = MPI_Wtime();

  if (list_compressed_blk[idx] != NULL){
    free(list_compressed_blk[idx]);
    list_compressed_blk[idx] = NULL;
  }
}

template <class Type>
void QubitRegister<Type>::DecompressTmpState()
{
  double t1, t2;
  int i;
  unsigned rank = openqu::mpi::Environment::rank();
  unsigned nprocs = openqu::mpi::Environment::size();
  MPI_Comm comm = openqu::mpi::Environment::comm();

  t1 = MPI_Wtime();

  if (tmp_state == NULL) {
    cerr << "unable to allocate memory for tmp_state, rank = " << rank << endl << flush;
  }
  BlkDecompress(tmp_compressed_blk, (double*) tmp_state, tmp_len);
  t2 = MPI_Wtime();
  if (tmp_compressed_blk != NULL){
    free(tmp_compressed_blk);
    tmp_compressed_blk = NULL;
  }
  
}

template <class Type>
int QubitRegister<Type>::CheckCache(int idx_1, unsigned char* b1, size_t len_1, int idx_2, unsigned char* b2, size_t len_2){

  if (enable_blk_cache == 0){
    return 2;
  }
  for(int i = 0; i < cache_top; i = i + 4){
    if(CmpBlk(b1, len_1, list_cache_blk[i], list_cache_len[i]) == 0 && CmpBlk(b2, len_2, list_cache_blk[i+1], list_cache_len[i+1]) == 0){
      CopyCompressedBlk(idx_1, i+2);
      CopyCompressedBlk(idx_2, i+3);
      blk_cache_hit++;
      return 0;
    }
  }
  blk_cache_miss++;

  if (cache_top < num_cache_line){
    list_cache_len[cache_top] = len_1;
    list_cache_blk[cache_top] = (unsigned char*) malloc(len_1);
    memcpy(list_cache_blk[cache_top], b1, len_1);
    list_cache_len[cache_top + 1] = len_2;
    list_cache_blk[cache_top + 1] = (unsigned char*) malloc(len_2);
    memcpy(list_cache_blk[cache_top + 1], b2, len_2);
    blk_cache_size = blk_cache_size + len_1 + len_2;
    return 1;
  }
  return 2;
}

template <class Type>
int QubitRegister<Type>::CheckCache(int idx_1, unsigned char* b1, size_t len_1){
  if (enable_blk_cache == 0){
    return 2;
  }
  for(int i = 0; i < cache_top; i = i + 4){
    if(CmpBlk(b1, len_1, list_cache_blk[i], list_cache_len[i]) == 0){
      CopyCompressedBlk(idx_1, i+2);
      blk_cache_hit++;
      return 0;
    }
  }
  blk_cache_miss++;

  if (cache_top < num_cache_line){
    list_cache_len[cache_top] = len_1;
    list_cache_blk[cache_top] = (unsigned char*) malloc(len_1);
    memcpy(list_cache_blk[cache_top], b1, len_1);
    blk_cache_size = blk_cache_size + len_1;
    return 1;
  }
  return 2;
}

template <class Type>
void QubitRegister<Type>::SetCache(int idx_1, unsigned char* b1, size_t len_1, int idx_2, unsigned char* b2, size_t len_2){
  list_cache_len[cache_top + 2] = len_1;
  list_cache_blk[cache_top + 2] = (unsigned char*) malloc(len_1);
  memcpy(list_cache_blk[cache_top + 2], b1, len_1);

  list_cache_len[cache_top + 3] = len_2;
  list_cache_blk[cache_top + 3] = (unsigned char*) malloc(len_2);
  memcpy(list_cache_blk[cache_top + 3], b2, len_2);

  cache_top = cache_top + 4;
  blk_cache_size = blk_cache_size + len_1 + len_2;
}

template <class Type>
void QubitRegister<Type>::SetCache(int idx_1, unsigned char* b1, size_t len_1){
  list_cache_len[cache_top + 2] = len_1;
  list_cache_blk[cache_top + 2] = (unsigned char*) malloc(len_1);
  memcpy(list_cache_blk[cache_top + 2], b1, len_1);

  cache_top = cache_top + 4;
  blk_cache_size = blk_cache_size + len_1;
}

template <class Type>
void QubitRegister<Type>::ClearCache(){
  for (int i = 0; i < num_cache_line * 4; i++){
    if (list_cache_blk[i] != NULL){
      free(list_cache_blk[i]);
      list_cache_blk[i] = NULL;
    } 
    list_cache_len[i] = 0;
  }
  cache_top = 0;
  if (hit_rate_0_count > 32 && (blk_cache_hit == 0)){
    enable_blk_cache = 0;
  }
}

template <class Type>
void QubitRegister<Type>::CopyCompressedBlk(int dest, int src){
  if (dest < 0){
    if (tmp_compressed_blk != NULL){
      free(tmp_compressed_blk);
      tmp_compressed_blk = NULL;
    }
    tmp_len = list_cache_len[src];
    tmp_compressed_blk = (unsigned char*) malloc(list_cache_len[src]);
    memcpy(tmp_compressed_blk, list_cache_blk[src], list_cache_len[src]);

  } else {
    if (list_compressed_blk[dest] != NULL){
      free(list_compressed_blk[dest]);
      list_compressed_blk[dest] = NULL;
    }
    list_len[dest] = list_cache_len[src];
    list_compressed_blk[dest] = (unsigned char*) malloc(list_cache_len[src]);
    memcpy(list_compressed_blk[dest], list_cache_blk[src], list_cache_len[src]);
  }
  
}

template <class Type>
int QubitRegister<Type>::CmpBlk(unsigned char* b1, size_t len_1, unsigned char* b2, size_t len_2){
  int res = 0;
  if (len_1 != len_2){
    res = 1;
  } else {
    for (int i = 0; i < len_1; i++){
      if (b1[i] != b2[i]){
        res = 1;
        break;
      }
    }
  }
  return res;

}

template <class Type>
void QubitRegister<Type>::OutputLclState(int idx)
{
  BlkDecompress(list_compressed_blk[idx], (double*) state, list_len[idx]);
}

template <class Type>
void QubitRegister<Type>::CompressionInfo()
{
  unsigned rank = openqu::mpi::Environment::rank();
  unsigned nprocs = openqu::mpi::Environment::size();
  MPI_Comm comm = openqu::mpi::Environment::comm();

  int mynoderank;
  MPI_Comm nodeComm;
  MPI_Comm_split_type( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank,
                       MPI_INFO_NULL, &nodeComm );
  MPI_Comm_rank(nodeComm, &mynoderank);

  double t;
  double sum_size = 0;
  double total_size = 0;
  double global_size = 0;
  double blk_cache_hit_rate = (blk_cache_hit + blk_cache_miss) ? blk_cache_hit / (blk_cache_hit + blk_cache_miss) : 0;
  double total_hit_rate;
  double max_compression_time, max_decompression_time;
  double max_compute_time, max_network_time;
  size_t total_blk_cache_size;

  for (int i = 0; i < num_block; i++){
    sum_size = sum_size + list_len[i];
  }

  t = sec();
  MPI_Reduce(&sum_size, &global_size, 1, MPI_DOUBLE , MPI_SUM, 0, comm);
  MPI_Reduce(&blk_cache_size, &total_blk_cache_size, 1, MPI_UINT64_T , MPI_SUM, 0, comm);

  MPI_Reduce(&blk_cache_hit_rate, &total_hit_rate, 1, MPI_DOUBLE , MPI_SUM, 0, comm);
  MPI_Reduce(&compression_time, &max_compression_time, 1, MPI_DOUBLE , MPI_MAX, 0, comm);
  MPI_Reduce(&decompression_time, &max_decompression_time, 1, MPI_DOUBLE , MPI_MAX, 0, comm);
  MPI_Reduce(&compute_time, &max_compute_time, 1, MPI_DOUBLE , MPI_MAX, 0, comm);
  network_time += sec() - t;
  MPI_Reduce(&network_time, &max_network_time, 1, MPI_DOUBLE , MPI_MAX, 0, comm);
  if (!rank) {
    cout << "Total compressed size (MB):\t" << global_size/1024/1024 << endl;
    cout << "Total block cache size (MB):\t" << total_blk_cache_size/1024/1024 << endl;
    cout << "Block cache hit rate:\t" << total_hit_rate  / nprocs << endl;
    cout << "Compression Time:\t" << compression_time << endl;
    cout << "Decompression Time:\t" << decompression_time << endl;
    cout << "Compute Time:\t" << compute_time << endl;
    cout << "Network Time:\t" << network_time << endl;
    cout << "Compress Level:\t" << compress_level << endl;
    
  }
  MPI_Reduce(&sum_size, &total_size, 1, MPI_DOUBLE , MPI_SUM, 0, nodeComm);
  MPI_Reduce(&blk_cache_size, &total_blk_cache_size, 1, MPI_UINT64_T , MPI_SUM, 0, nodeComm);
  if (!mynoderank){
    //cout << "Node " << rank / x_rpn << ": Total compressed size (MB):\t" << total_size/1024/1024 << endl;
    //cout << "Node " << rank / x_rpn << ": Total block cache size (MB):\t" << total_blk_cache_size/1024/1024 << endl;
    if (compress_level < 2) {
      if ((total_size+total_blk_cache_size)/1024/1024 > 32*1024){
        compress_level = compress_level + 1;
      }
    }
    else if (compress_level < 3) {
      if ((total_size+total_blk_cache_size)/1024/1024 > 64*1024){
        compress_level = compress_level + 1;
      }
    }
    else if (compress_level < 4) {
      if ((total_size+total_blk_cache_size)/1024/1024 > 96*1024){
        compress_level = compress_level + 1;
      }
    }
    else if (compress_level < LEVEL - 1) {
      if ((total_size+total_blk_cache_size)/1024/1024 > 128*1024){
        compress_level = compress_level + 1;
      }
    }
  }
  if (enable_blk_cache && (blk_cache_hit_rate == 0) && (blk_cache_miss > 0)){
    hit_rate_0_count++;
  }

  compression_time = 0;
  decompression_time = 0;
  compute_time = 0;
  network_time = 0;
  blk_cache_size = 0;
  blk_cache_hit = 0;
  blk_cache_miss = 0;
  MPI_Bcast(&compress_level, 1 , MPI_INT, 0, nodeComm);
}

template <class Type>
void QubitRegister<Type>::OutputCompressedByteToFile()
{
  unsigned rank = openqu::mpi::Environment::rank();
  unsigned nprocs = openqu::mpi::Environment::size();
  MPI_Comm comm = openqu::mpi::Environment::comm();

  char filename[256];
  size_t sum_size = 0;
  size_t total_size = 0;

  double t0 = sec();
  openqu::mpi::barrier();
  for (int i = 0; i < num_block; i++){
    sprintf(filename, "szcomp-%u-%d", rank, i);
    FILE *pFile = fopen(filename, "wb");
    if (pFile == NULL){
      cerr << "Failed to open the file: " << filename << endl;
    } else {
    fwrite(list_compressed_blk[i], 1, list_len[i], pFile);
    fclose(pFile);
    sum_size = sum_size + list_len[i];
    }
  }
  sprintf(filename, "compress_level-%u", rank);
  FILE *pFile = fopen(filename, "w");
  if (pFile == NULL){
    cerr << "Failed to open the file: " << filename << endl;
  } else {
    fprintf(pFile, "%d\n", compress_level);
    fclose(pFile);
  }

  openqu::mpi::barrier();
  double t1 = sec();

  MPI_Reduce(&sum_size, &total_size, 1, MPI_UINT64_T , MPI_SUM, 0, comm);
  if (!rank) {
    cout << "Total output compressed size (MB):\t" << total_size/1024/1024 << endl;
    cout << "Output CB time:\t" << t1 - t0 << endl;
  }

}

template <class Type>
void QubitRegister<Type>::InputCompressedByteFromFile()
{
  unsigned rank = openqu::mpi::Environment::rank();
  unsigned nprocs = openqu::mpi::Environment::size();
  MPI_Comm comm = openqu::mpi::Environment::comm();

  char filename[256];
  size_t sum_size = 0;
  size_t total_size = 0;

  double t0 = sec();
  openqu::mpi::barrier();

  for (int i = 0; i < num_block; i++){
    sprintf(filename, "szcomp-%u-%d", rank, i);
    FILE *pFile = fopen(filename, "rb");
    if (pFile == NULL){
      cerr << "No Input file: " << filename << endl;
    } else {
      fseek(pFile, 0, SEEK_END);
      list_len[i] = ftell(pFile);
      if (list_compressed_blk[i] != NULL) free(list_compressed_blk[i]);
      list_compressed_blk[i] = (unsigned char*) malloc(list_len[i]);
      fseek(pFile, 0, SEEK_SET);
      fread(list_compressed_blk[i], 1, list_len[i], pFile);
      fclose(pFile);
      sum_size = sum_size + list_len[i];
    }
  }

  sprintf(filename, "compress_level-%u", rank);
  FILE *pFile = fopen(filename, "r");
  if (pFile == NULL){
    cerr << "Failed to open the file: " << filename << "set compress_level to 2"<< endl;
    compress_level = 2;
  } else {
    fscanf(pFile, "%d", &compress_level);
    fclose(pFile);
  }

  openqu::mpi::barrier();
  double t1 = sec();

  MPI_Reduce(&sum_size, &total_size, 1, MPI_UINT64_T , MPI_SUM, 0, comm);
  if (!rank) {
    cout << "Total input compressed size (MB):\t" << total_size/1024/1024 << endl;
    cout << "Input CB time:\t" << t1 - t0 << endl;
  }
}


