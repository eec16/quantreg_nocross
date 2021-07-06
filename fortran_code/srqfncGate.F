C==============================================================
C     srqfncGate.F
C
C     Gateway for srqfnc.f from no cross quant reg algorithm
C     By Evan Corden 7/2/2021
C     This is a MEX-file for MATLAB.
C==============================================================
#include "fintrf.h"
C     The gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

C--------------------------------------------------------------
C Declare initial arguments and function return types
      implicit none
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
      double precision tmp_int 
      double precision, allocatable, dimension(:) :: inta_tmp 
      double precision, allocatable, dimension(:,:) :: intm_tmp 
      mwPointer mxCreateDoubleMatrix, mxGetPr
      mwPointer mxGetM, mxGetN
      mwPointer mdim, n, sol_n
      mwSize size
      
C--------------------------------------------------------------
C Single integer values
C nnza1
        mwPointer nnza1_pr
        integer nnza1
C nnza2
        mwPointer nnza2_pr
        integer nnza2
C m
        mwPointer m_pr
        integer m
C n1
        mwPointer n1_pr
        integer n1
C n2
        mwPointer n2_pr
        integer n2
C nnzdmax
        mwPointer nnzdmax_pr
        integer nnzdmax
C nnzemax
        mwPointer nnzemax_pr
        integer nnzemax
C nnzgmax
        mwPointer nnzgmax_pr
        integer nnzgmax
C nnzhmax
        mwPointer nnzhmax_pr
        integer nnzhmax
C iwmax
        mwPointer iwmax_pr
        integer iwmax
C nnzlmax
        mwPointer nnzlmax_pr
        integer nnzlmax
C nsubmax
        mwPointer nsubmax_pr
        integer nsubmax
C cachsz
        mwPointer cachsz_pr
        integer cachsz
C level
        mwPointer level_pr
        integer level
C tmpmax
        mwPointer tmpmax_pr
        integer tmpmax
C ierr
        mwPointer ierr_pr
        integer ierr
C maxit
        mwPointer maxit_pr
        integer maxit
C maxn1n2
        mwPointer maxn1n2_pr
        integer maxn1n2
        
C-------------------------------------------------------------------------
C Integer arrays
C-------------------------------------------------------------------------
C ja1
        mwPointer ja1_pr
        integer, allocatable, dimension(:) :: ja1
C jao1
        mwPointer jao1_pr
        integer, allocatable, dimension(:) :: jao1
C ja2
        mwPointer ja2_pr
        integer, allocatable, dimension(:) :: ja2
C jao2
        mwPointer jao2_pr
        integer, allocatable, dimension(:) :: jao2
C jdsub
        mwPointer jdsub_pr
        integer, allocatable, dimension(:) :: jdsub
C jd
        mwPointer jd_pr
        integer, allocatable, dimension(:) :: jd
C ia1
        mwPointer ia1_pr
        integer, allocatable, dimension(:) :: ia1
C iao1
        mwPointer iao1_pr
        integer, allocatable, dimension(:) :: iao1
C ia2
        mwPointer ia2_pr
        integer, allocatable, dimension(:) :: ia2
C iao2
        mwPointer iao2_pr
        integer, allocatable, dimension(:) :: iao2
C id
        mwPointer id_pr
        integer, allocatable, dimension(:) :: id
C lindx
        mwPointer lindx_pr
        integer, allocatable, dimension(:) :: lindx
C xlindx
        mwPointer xlindx_pr
        integer, allocatable, dimension(:) :: xlindx
C iw
        mwPointer iw_pr
        integer, allocatable, dimension(:,:) :: iw
C xlnz
        mwPointer xlnz_pr
        integer, allocatable, dimension(:) :: xlnz
C iwork
        mwPointer iwork_pr
        integer, allocatable, dimension(:) :: iwork
C xsuper
        mwPointer xsuper_pr
        integer, allocatable, dimension(:) :: xsuper
C je
        mwPointer je_pr
        integer, allocatable, dimension(:) :: je
C ie
        mwPointer ie_pr
        integer, allocatable, dimension(:) :: ie
C jg
        mwPointer jg_pr
        integer, allocatable, dimension(:) :: jg
C ig
        mwPointer ig_pr
        integer, allocatable, dimension(:) :: ig
C jh
        mwPointer jh_pr
        integer, allocatable, dimension(:) :: jh
C ih
        mwPointer ih_pr
        integer, allocatable, dimension(:) :: ih
C-------------------------------------------------------------------------
C Doubles
C-------------------------------------------------------------------------
C small
        mwPointer small_pr
        double precision small
C a1
        mwPointer a1_pr
        double precision, allocatable, dimension(:) :: a1
C ao1
        mwPointer ao1_pr
        double precision, allocatable, dimension(:) :: ao1
C a2
        mwPointer a2_pr
        double precision, allocatable, dimension(:) :: a2
C ao2
        mwPointer ao2_pr
        double precision, allocatable, dimension(:) :: ao2
C dsub
        mwPointer dsub_pr
        double precision, allocatable, dimension(:) :: dsub
C d
        mwPointer d_pr
        double precision, allocatable, dimension(:) :: d
C g
        mwPointer g_pr
        double precision, allocatable, dimension(:) :: g
C h
        mwPointer h_pr
        double precision, allocatable, dimension(:) :: h
C lnz
        mwPointer lnz_pr
        double precision, allocatable, dimension(:) :: lnz
C c1
        mwPointer c1_pr
        double precision, allocatable, dimension(:) :: c1
C c2
        mwPointer c2_pr
        double precision, allocatable, dimension(:) :: c2
C y
        mwPointer y_pr
        double precision, allocatable, dimension(:) :: y
C ww1
        mwPointer ww1_pr
        double precision, allocatable, dimension(:) :: ww1
C wwm
        mwPointer wwm_pr
        double precision, allocatable, dimension(:,:) :: wwm
C tmpvec
        mwPointer tmpvec_pr
        double precision, allocatable, dimension(:) :: tmpvec
C wwn1
        mwPointer wwn1_pr
        double precision, allocatable, dimension(:,:) :: wwn1
C wwn2
        mwPointer wwn2_pr
        double precision, allocatable, dimension(:,:) :: wwn2
C x1
        mwPointer x1_pr
        double precision, allocatable, dimension(:) :: x1
C x2
        mwPointer x2_pr
        double precision, allocatable, dimension(:) :: x2
C s
        mwPointer s_pr
        double precision, allocatable, dimension(:) :: s
C u
        mwPointer u_pr
        double precision, allocatable, dimension(:) :: u
C e
        mwPointer e_pr
        double precision, allocatable, dimension(:) :: e
C timewd
        mwPointer timewd_pr
        double precision, allocatable, dimension(:) :: timewd
C-------------------------------------------------------------------------
C Output declarations
C-------------------------------------------------------------------------
        mwPointer y_coeff_pr, y_maxiter_pr, y_ierr_pr, y_timewd_pr
        mwPointer test_pr
        integer i

C-------------------------------------------------------------------------
C Statements
C-------------------------------------------------------------------------
C Single Integers
C-------------------------------------------------------------------------
C nnza1
        nnza1_pr = mxGetPr(prhs(3))
        call mxCopyPtrToReal8(nnza1_pr, tmp_int, 1)
        nnza1 = IDINT(tmp_int)
        
C nnza2
        nnza2_pr = mxGetPr(prhs(11))
        call mxCopyPtrToReal8(nnza2_pr, tmp_int, 1)
        nnza2 = IDINT(tmp_int)
C m
        m_pr = mxGetPr(prhs(2))
        call mxCopyPtrToReal8(m_pr, tmp_int, 1)
        m = IDINT(tmp_int)
C n1
        n1_pr = mxGetPr(prhs(1))
        call mxCopyPtrToReal8(n1_pr, tmp_int, 1)
        n1 = IDINT(tmp_int)
C n2
        n2_pr = mxGetPr(prhs(10))
        call mxCopyPtrToReal8(n2_pr, tmp_int, 1)
        n2 = IDINT(tmp_int)
C nnzdmax
        nnzdmax_pr = mxGetPr(prhs(18))
        call mxCopyPtrToReal8(nnzdmax_pr, tmp_int, 1)
        nnzdmax = IDINT(tmp_int)
C nnzemax
        nnzemax_pr = mxGetPr(prhs(24))
        call mxCopyPtrToReal8(nnzemax_pr, tmp_int, 1)
        nnzemax = IDINT(tmp_int)
C nnzgmax
        nnzgmax_pr = mxGetPr(prhs(28))
        call mxCopyPtrToReal8(nnzgmax_pr, tmp_int, 1)
        nnzgmax = IDINT(tmp_int)
C nnzhmax
        nnzhmax_pr = mxGetPr(prhs(32))
        call mxCopyPtrToReal8(nnzhmax_pr, tmp_int, 1)
        nnzhmax = IDINT(tmp_int)
C iwmax
        iwmax_pr = mxGetPr(prhs(43))
        call mxCopyPtrToReal8(iwmax_pr, tmp_int, 1)
        iwmax = IDINT(tmp_int)
C nnzlmax
        nnzlmax_pr = mxGetPr(prhs(39))
        call mxCopyPtrToReal8(nnzlmax_pr, tmp_int, 1)
        nnzlmax = IDINT(tmp_int)
C nsubmax
        nsubmax_pr = mxGetPr(prhs(36))
        call mxCopyPtrToReal8(nsubmax_pr, tmp_int, 1)
        nsubmax = IDINT(tmp_int)
C cachsz
        cachsz_pr = mxGetPr(prhs(53))
        call mxCopyPtrToReal8(cachsz_pr, tmp_int, 1)
        cachsz = IDINT(tmp_int)
C level
        level_pr = mxGetPr(prhs(54))
        call mxCopyPtrToReal8(level_pr, tmp_int, 1)
        level = IDINT(tmp_int)
C tmpmax
        tmpmax_pr = mxGetPr(prhs(46))
        call mxCopyPtrToReal8(tmpmax_pr, tmp_int, 1)
        tmpmax = IDINT(tmp_int)
C ierr
        ierr_pr = mxGetPr(prhs(63))
        call mxCopyPtrToReal8(ierr_pr, tmp_int, 1)
        ierr = IDINT(tmp_int)
C maxit
        maxit_pr = mxGetPr(prhs(64))
        call mxCopyPtrToReal8(maxit_pr, tmp_int, 1)
        maxit = IDINT(tmp_int)
C maxn1n2
        maxn1n2_pr = mxGetPr(prhs(48))
        call mxCopyPtrToReal8(maxn1n2_pr, tmp_int, 1)
        maxn1n2 = IDINT(tmp_int)
C-------------------------------------------------------------------------
C Integer arrays
C-------------------------------------------------------------------------
C ja1
        n = mxGetN(prhs(5))
        ja1_pr = mxGetPr(prhs(5))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(ja1_pr, inta_tmp, n)
        allocate(ja1(n))
        ja1 = IDINT(inta_tmp)
        deallocate(inta_tmp)
C jao1
        n = mxGetN(prhs(8))
        jao1_pr = mxGetPr(prhs(8))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(jao1_pr, inta_tmp, n)
        allocate(jao1(n))
        jao1 = IDINT(inta_tmp)
        deallocate(inta_tmp)
C ja2
        n = mxGetN(prhs(13))
        ja2_pr = mxGetPr(prhs(13))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(ja2_pr, inta_tmp, n)
        allocate(ja2(n))
        ja2 = IDINT(inta_tmp)
        deallocate(inta_tmp)
C jao2
        n = mxGetN(prhs(16))
        jao2_pr = mxGetPr(prhs(16))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(jao2_pr, inta_tmp, n)
        allocate(jao2(n))
        jao2 = IDINT(inta_tmp)
        deallocate(inta_tmp)
C jdsub
        n = mxGetN(prhs(23))
        jdsub_pr = mxGetPr(prhs(23))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(jdsub_pr, inta_tmp, n)
        allocate(jdsub(n))
        jdsub = IDINT(inta_tmp)
        deallocate(inta_tmp)
C jd
        n = mxGetN(prhs(20))
        jd_pr = mxGetPr(prhs(20))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(jd_pr, inta_tmp, n)
        allocate(jd(n))
        jd = IDINT(inta_tmp)
        deallocate(inta_tmp)
C ia1
        n = mxGetN(prhs(6))
        ia1_pr = mxGetPr(prhs(6))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(ia1_pr, inta_tmp, n)
        allocate(ia1(n))
        ia1 = IDINT(inta_tmp)
        deallocate(inta_tmp)
C iao1
        n = mxGetN(prhs(9))
        iao1_pr = mxGetPr(prhs(9))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(iao1_pr, inta_tmp, n)
        allocate(iao1(n))
        iao1 = IDINT(inta_tmp)
        deallocate(inta_tmp)
C ia2
        n = mxGetN(prhs(14))
        ia2_pr = mxGetPr(prhs(14))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(ia2_pr, inta_tmp, n)
        allocate(ia2(n))
        ia2 = IDINT(inta_tmp)
        deallocate(inta_tmp)
C iao2
        n = mxGetN(prhs(17))
        iao2_pr = mxGetPr(prhs(17))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(iao2_pr, inta_tmp, n)
        allocate(iao2(n))
        iao2 = IDINT(inta_tmp)
        deallocate(inta_tmp)
C id
        n = mxGetN(prhs(21))
        id_pr = mxGetPr(prhs(21))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(id_pr, inta_tmp, n)
        allocate(id(n))
        id = IDINT(inta_tmp)
        deallocate(inta_tmp)
C lindx
        n = mxGetN(prhs(37))
        lindx_pr = mxGetPr(prhs(37))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(lindx_pr, inta_tmp, n)
        allocate(lindx(n))
        lindx = IDINT(inta_tmp)
        deallocate(inta_tmp)
C xlindx
        n = mxGetN(prhs(38))
        xlindx_pr = mxGetPr(prhs(38))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(xlindx_pr, inta_tmp, n)
        allocate(xlindx(n))
        xlindx = IDINT(inta_tmp)
        deallocate(inta_tmp)
C iw
        n = mxGetN(prhs(42))
        mdim = mxGetM(prhs(42))
        size = mdim*n
        iw_pr = mxGetPr(prhs(42))
        allocate(intm_tmp(mdim,n))
        call mxCopyPtrToReal8(iw_pr, intm_tmp, size)
        allocate(iw(mdim,n))
        iw = IDINT(intm_tmp)
        deallocate(intm_tmp)
C xlnz
        n = mxGetN(prhs(41))
        xlnz_pr = mxGetPr(prhs(41))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(xlnz_pr, inta_tmp, n)
        allocate(xlnz(n))
        xlnz = IDINT(inta_tmp)
        deallocate(inta_tmp)
C iwork
        n = mxGetN(prhs(44))
        iwork_pr = mxGetPr(prhs(44))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(iwork_pr, inta_tmp, n)
        allocate(iwork(n))
        iwork = IDINT(inta_tmp)
        deallocate(inta_tmp)
C xsuper
        n = mxGetN(prhs(45))
        xsuper_pr = mxGetPr(prhs(45))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(xsuper_pr, inta_tmp, n)
        allocate(xsuper(n))
        xsuper = IDINT(inta_tmp)
        deallocate(inta_tmp)
C je
        n = mxGetN(prhs(26))
        je_pr = mxGetPr(prhs(26))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(je_pr, inta_tmp, n)
        allocate(je(n))
        je = IDINT(inta_tmp)
        deallocate(inta_tmp)
C ie
        n = mxGetN(prhs(27))
        ie_pr = mxGetPr(prhs(27))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(ie_pr, inta_tmp, n)
        allocate(ie(n))
        ie = IDINT(inta_tmp)
        deallocate(inta_tmp)
C jg
        n = mxGetN(prhs(30))
        jg_pr = mxGetPr(prhs(30))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(jg_pr, inta_tmp, n)
        allocate(jg(n))
        jg = IDINT(inta_tmp)
        deallocate(inta_tmp)
C ig
        n = mxGetN(prhs(31))
        ig_pr = mxGetPr(prhs(31))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(ig_pr, inta_tmp, n)
        allocate(ig(n))
        ig = IDINT(inta_tmp)
        deallocate(inta_tmp)
C jh
        n = mxGetN(prhs(34))
        jh_pr = mxGetPr(prhs(34))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(jh_pr, inta_tmp, n)
        allocate(jh(n))
        jh = IDINT(inta_tmp)
        deallocate(inta_tmp)
C ih
        n = mxGetN(prhs(35))
        ih_pr = mxGetPr(prhs(35))
        allocate(inta_tmp(n))
        call mxCopyPtrToReal8(ih_pr, inta_tmp, n)
        allocate(ih(n))
        ih = IDINT(inta_tmp)
        deallocate(inta_tmp)
C-------------------------------------------------------------------------
C Doubles
C-------------------------------------------------------------------------      
C small
        small_pr = mxGetPr(prhs(62))
        call mxCopyPtrToReal8(small_pr, small, 1)
C a1
        n = mxGetN(prhs(4))
        a1_pr = mxGetPr(prhs(4))
        allocate(a1(n))
        call mxCopyPtrToReal8(a1_pr, a1, n)
C ao1
        n = mxGetN(prhs(7))
        ao1_pr = mxGetPr(prhs(7))
        allocate(ao1(n))
        call mxCopyPtrToReal8(ao1_pr, ao1, n)
C a2
        n = mxGetN(prhs(12))
        a2_pr = mxGetPr(prhs(12))
        allocate(a2(n))
        call mxCopyPtrToReal8(a2_pr, a2, n)
C ao2
        n = mxGetN(prhs(15))
        ao2_pr = mxGetPr(prhs(15))
        allocate(ao2(n))
        call mxCopyPtrToReal8(ao2_pr, ao2, n)
C dsub
        n = mxGetN(prhs(22))
        dsub_pr = mxGetPr(prhs(22))
        allocate(dsub(n))
        call mxCopyPtrToReal8(dsub_pr, dsub, n)
C d
        n = mxGetN(prhs(19))
        d_pr = mxGetPr(prhs(19))
        allocate(d(n))
        call mxCopyPtrToReal8(d_pr, d, n)
C g
        n = mxGetN(prhs(29))
        g_pr = mxGetPr(prhs(29))
        allocate(g(n))
        call mxCopyPtrToReal8(g_pr, g, n)
C h
        n = mxGetN(prhs(33))
        h_pr = mxGetPr(prhs(33))
        allocate(h(n))
        call mxCopyPtrToReal8(h_pr, h, n)
C lnz
        n = mxGetN(prhs(40))
        lnz_pr = mxGetPr(prhs(40))
        allocate(lnz(n))
        call mxCopyPtrToReal8(lnz_pr, lnz, n)
C c1
        n = mxGetN(prhs(59))
        c1_pr = mxGetPr(prhs(59))
        allocate(c1(n))
        call mxCopyPtrToReal8(c1_pr, c1, n)
C c2
        n = mxGetN(prhs(60))
        c2_pr = mxGetPr(prhs(60))
        allocate(c2(n))
        call mxCopyPtrToReal8(c2_pr, c2, n)
C y
        sol_n = mxGetN(prhs(61))
        y_pr = mxGetPr(prhs(61))
        allocate(y(sol_n))
        call mxCopyPtrToReal8(y_pr, y, sol_n)
C ww1
        n = mxGetN(prhs(49))
        ww1_pr = mxGetPr(prhs(49))
        allocate(ww1(n))
        call mxCopyPtrToReal8(ww1_pr, ww1, n)
C wwm
        mdim = mxGetM(prhs(50))
        n = mxGetN(prhs(50))
        size = mdim*n
        wwm_pr = mxGetPr(prhs(50))
        allocate(wwm(mdim,n))
        call mxCopyPtrToReal8(wwm_pr, wwm, size)
C tmpvec
        n = mxGetN(prhs(47))
        tmpvec_pr = mxGetPr(prhs(47))
        allocate(tmpvec(n))
        call mxCopyPtrToReal8(tmpvec_pr, tmpvec, n)
C wwn1
        mdim = mxGetM(prhs(51))
        n = mxGetN(prhs(51))
        size = mdim*n
        wwn1_pr = mxGetPr(prhs(51))
        allocate(wwn1(mdim,n))
        call mxCopyPtrToReal8(wwn1_pr, wwn1, size)

C wwn2
        mdim = mxGetM(prhs(52))
        n = mxGetN(prhs(52))
        size = mdim*n
        wwn2_pr = mxGetPr(prhs(52))
        allocate(wwn2(mdim,n))
        call mxCopyPtrToReal8(wwn2_pr, wwn2, size)
C x1
        n = mxGetN(prhs(55))
        x1_pr = mxGetPr(prhs(55))
        allocate(x1(n))
        call mxCopyPtrToReal8(x1_pr, x1, n)
C x2
        n = mxGetN(prhs(56))
        x2_pr = mxGetPr(prhs(56))
        allocate(x2(n))
        call mxCopyPtrToReal8(x2_pr, x2, n)
C s
        n = mxGetN(prhs(57))
        s_pr = mxGetPr(prhs(57))
        allocate(s(n))
        call mxCopyPtrToReal8(s_pr, s, n)
C u
        n = mxGetN(prhs(58))
        u_pr = mxGetPr(prhs(58))
        allocate(u(n))
        call mxCopyPtrToReal8(u_pr, u, n)
C e
        n = mxGetN(prhs(25))
        e_pr = mxGetPr(prhs(25))
        allocate(e(n))
        call mxCopyPtrToReal8(e_pr, e, n)
C timewd
        n = mxGetN(prhs(65))
        timewd_pr = mxGetPr(prhs(65))
        allocate(timewd(n))
        call mxCopyPtrToReal8(timewd_pr, timewd, n) 

C     Call the computational subroutine
        call srqfnc(n1,m,nnza1,a1,ja1,ia1,ao1,jao1,iao1,n2,nnza2,
     &                  a2,ja2,ia2,ao2,jao2,iao2,nnzdmax,d,jd,id,
     &                  dsub,jdsub,nnzemax,e,je,ie,nnzgmax,g,jg,ig,
     &                  nnzhmax,h,jh,ih,nsubmax,lindx,xlindx,nnzlmax,
     &                  lnz,xlnz,iw,iwmax,iwork,xsuper,tmpmax,tmpvec,
     &                  maxn1n2,ww1,wwm,wwn1,wwn2,cachsz,level,x1,x2,
     &                  s,u,c1,c2,y,small,ierr,maxit,timewd)

C     Deallocate dynamic memory
        deallocate(ja1)
        deallocate(jao1)
        deallocate(ja2)
        deallocate(jao2)
        deallocate(jdsub)
        deallocate(jd)
        deallocate(ia1)
        deallocate(iao1)
        deallocate(ia2)
        deallocate(iao2)
        deallocate(id)
        deallocate(lindx)
        deallocate(xlindx)
        deallocate(iw)
        deallocate(xlnz)
        deallocate(iwork)
        deallocate(xsuper)
        deallocate(je)
        deallocate(ie)
        deallocate(jg)
        deallocate(ig)
        deallocate(jh)
        deallocate(ih)
        deallocate(a1)
        deallocate(ao1)
        deallocate(a2)
        deallocate(ao2)
        deallocate(dsub)
        deallocate(d)
        deallocate(g)
        deallocate(h)
        deallocate(lnz)
        deallocate(c1)
        deallocate(c2)
        deallocate(ww1)
        deallocate(wwm)
        deallocate(tmpvec)
        deallocate(wwn1)
        deallocate(wwn2)
        deallocate(x1)
        deallocate(x2)
        deallocate(s)
        deallocate(u)
        deallocate(e)
          
C     Load the output data
        plhs(1) = mxCreateDoubleMatrix(1,sol_n,0)
        plhs(2) = mxCreateDoubleMatrix(1,1,0)
        plhs(3) = mxCreateDoubleMatrix(1,1,0)
        plhs(4) = mxCreateDoubleMatrix(1,7,0)

        y_coeff_pr = mxGetPr(plhs(1))
        y_ierr_pr = mxGetPr(plhs(2))
        y_maxiter_pr = mxGetPr(plhs(3))
        y_timewd_pr = mxGetPr(plhs(4))

        call mxCopyReal8ToPtr(y,y_coeff_pr,sol_n)
        tmp_int = DBLE(ierr)
        call mxCopyReal8ToPtr(tmp_int,y_ierr_pr,1)
        tmp_int = DBLE(maxit)
        call mxCopyReal8ToPtr(tmp_int,y_maxiter_pr,1)
        call mxCopyReal8ToPtr(timewd,y_timewd_pr,7)

C Deallocate 2 dynamic outputs
        deallocate(y)
        deallocate(timewd)
      return
      end
