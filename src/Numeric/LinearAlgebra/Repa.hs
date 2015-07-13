{-# LANGUAGE FlexibleContexts #-}

module Numeric.LinearAlgebra.Repa 
  ( Numeric
  , Field
  , Product
  , HShape(..)
  , LSDiv
  -- * Dot product
  , dot
  , dotS
  , dotSIO
  , dotP
  , dotPIO
  -- * Dense matrix-vector product.
  , app
  , appS
  , appSIO
  , appP
  , appPIO
  -- * Dense matrix-matrix product.
  , mul
  , mulS
  , mulSIO
  , mulP
  , mulPIO
  -- * Vector outer product.
  , outer
  , outerS
  , outerSIO
  , outerP
  , outerPIO
  -- * Kronecker product.
  , kronecker
  , kroneckerS
  , kroneckerSIO
  , kroneckerP
  , kroneckerPIO
  -- * Cross product.
  , cross
  , crossS
  , crossSIO
  -- * Sum of elements.
  , sumElements
  , sumElementsS
  , sumElementsSIO
  , sumElementsP
  , sumElementsPIO
  -- * Product of elements.
  , prodElements
  , prodElementsS
  , prodElementsSIO
  , prodElementsP
  , prodElementsPIO
  -- * Linear systems.
  , (<\>)
  , solve
  , solveS
  , solveSIO
  , solveP
  , solvePIO
  , linearSolve
  , linearSolveS
  , linearSolveSIO
  , linearSolveP
  , linearSolvePIO
  , linearSolveLS
  , linearSolveLS_S
  , linearSolveLS_SIO
  , linearSolveLS_P
  , linearSolveLS_PIO
  , linearSolveSVD
  , linearSolveSVD_S
  , linearSolveSVD_SIO
  , linearSolveSVD_P
  , linearSolveSVD_PIO
  , luSolve
  , luSolveS
  , luSolveSIO
  , luSolveP
  , luSolvePIO
  , cholSolve
  , cholSolveS
  , cholSolveSIO
  , cholSolveP
  , cholSolvePIO
  -- * Inverse and pseudoinverse
  , inv
  , invS
  , invSIO
  , invP
  , invPIO
  , pinv
  , pinvS
  , pinvSIO
  , pinvP
  , pinvPIO
  , pinvTol
  , pinvTolS
  , pinvTolSIO
  , pinvTolP
  , pinvTolPIO
  -- * Determinant and rank
  , rcond
  , rcondS
  , rcondSIO
  , rcondP
  , rcondPIO
  , rank
  , rankS
  , rankSIO
  , rankP
  , rankPIO
  , det
  , detS
  , detSIO
  , detP
  , detPIO
  , invlndet
  , invlndetS
  , invlndetSIO
  , invlndetP
  , invlndetPIO
  -- * Norms
  , norm_Frob
  , norm_FrobS
  , norm_FrobSIO
  , norm_FrobP
  , norm_FrobPIO
  , norm_nuclear
  , norm_nuclearS
  , norm_nuclearSIO
  , norm_nuclearP
  , norm_nuclearPIO
  -- * Nullspace and range
  , orth
  , orthS
  , orthSIO
  , orthP
  , orthPIO
  , nullspace
  , nullspaceS
  , nullspaceSIO
  , nullspaceP
  , nullspacePIO
  , null1
  , null1S
  , null1SIO
  , null1P
  , null1PIO
  , null1sym
  , null1symS
  , null1symSIO
  , null1symP
  , null1symPIO
  -- * SVD
  , svd
  , svdS
  , svdSIO
  , svdP
  , svdPIO
  , thinSVD
  , thinSVD_S
  , thinSVD_SIO
  , thinSVD_P
  , thinSVD_PIO
  , compactSVD
  , compactSVD_S
  , compactSVD_SIO
  , compactSVD_P
  , compactSVD_PIO
  , singularValues
  , singularValuesS
  , singularValuesSIO
  , singularValuesP
  , singularValuesPIO
  , leftSV
  , leftSV_S
  , leftSV_SIO
  , leftSV_P
  , leftSV_PIO
  , rightSV
  , rightSV_S
  , rightSV_SIO
  , rightSV_P
  , rightSV_PIO
  -- * Eigensystems
  , eig
  , eigS
  , eigSIO
  , eigP
  , eigPIO
  , eigSH
  , eigSH_S
  , eigSH_SIO
  , eigSH_P
  , eigSH_PIO
  , eigSH'
  , eigSH'S
  , eigSH'SIO
  , eigSH'P
  , eigSH'PIO
  , eigenvalues
  , eigenvaluesS
  , eigenvaluesSIO
  , eigenvaluesP
  , eigenvaluesPIO
  , eigenvaluesSH
  , eigenvaluesSH_S
  , eigenvaluesSH_SIO
  , eigenvaluesSH_P
  , eigenvaluesSH_PIO
  , eigenvaluesSH'
  , eigenvaluesSH'S
  , eigenvaluesSH'SIO
  , eigenvaluesSH'P
  , eigenvaluesSH'PIO
  , geigSH'
  , geigSH'S
  , geigSH'SIO
  , geigSH'P
  , geigSH'PIO
  -- * QR
  , qr
  , qrS
  , qrSIO
  , qrP
  , qrPIO
  , rq
  , rqS
  , rqSIO
  , rqP
  , rqPIO
  , qrRaw
  , qrRawS
  , qrRawSIO
  , qrRawP
  , qrRawPIO
  , qrgr
  -- * Cholesky
  , chol
  , cholS
  , cholSIO
  , cholP
  , cholPIO
  , chol'
  , chol'S
  , chol'SIO
  , chol'P
  , chol'PIO
  -- * Hessenberg
  , hess
  , hessS
  , hessSIO
  , hessP
  , hessPIO
  -- * Schur
  , schur
  , schurS
  , schurSIO
  , schurP
  , schurPIO
  ) where

import Numeric.LinearAlgebra.Repa.Conversion

import Data.Array.Repa hiding (rank)
import Data.Array.Repa.Repr.ForeignPtr
import qualified Numeric.LinearAlgebra.HMatrix as H
import Numeric.LinearAlgebra.HMatrix (Complex, Numeric, Field, LSDiv, Normed, Product, Vector)

-- Dot product

dot :: Numeric t => Array F DIM1 t -> Array F DIM1 t -> t
-- ^Vector dot product.
dot v u = repa2hv v `H.dot` repa2hv u 

dotS :: Numeric t => Array D DIM1 t -> Array D DIM1 t -> t
-- ^Vector dot product. Arguments computed sequentially.
dotS v u = repa2hvS v `H.dot` repa2hvS u

dotSIO :: Numeric t => Array D DIM1 t -> Array D DIM1 t -> IO t
-- ^Vector dot product. Arguments computed sequentially inside the IO monad.
dotSIO v u = H.dot <$> repa2hvSIO v <*> repa2hvSIO u

dotP :: (Numeric t, Monad m) => Array D DIM1 t -> Array D DIM1 t -> m t
-- ^Vector dot product. Arguments computed in parallel.
dotP v u = H.dot <$> repa2hvP v <*> repa2hvP u

dotPIO :: Numeric t => Array D DIM1 t -> Array D DIM1 t -> IO t
-- ^Vector dot product. Arguments computed in parallel inside the IO monad.
dotPIO v u = H.dot <$> repa2hvPIO v <*> repa2hvPIO u

-- Dense matrix-vector product

app :: Numeric t => Array F DIM2 t -> Array F DIM1 t -> Array F DIM1 t
-- ^Dense matrix-vector product.
app m v = hv2repa $ repa2hm m `H.app` repa2hv v

appS :: Numeric t => Array D DIM2 t -> Array D DIM1 t -> Array F DIM1 t
-- ^Dense matrix-vector product. Arguments computed sequentially.
appS m v = hv2repa $ repa2hmS m `H.app` repa2hvS v

appSIO :: Numeric t => Array D DIM2 t -> Array D DIM1 t -> IO (Array F DIM1 t)
-- ^Dense matrix-vector product. Arguments computed sequentially inside the IO monad.
appSIO m v = hv2repa <$> (H.app <$> repa2hmSIO m <*> repa2hvSIO v)

appP :: (Numeric t, Monad m) => Array D DIM2 t -> Array D DIM1 t -> m (Array F DIM1 t)
-- ^Dense matrix-vector product. Arguments computed in parallel.
appP m v = hv2repa <$> (H.app <$> repa2hmP m <*> repa2hvP v)

appPIO :: Numeric t => Array D DIM2 t -> Array D DIM1 t -> IO (Array F DIM1 t)
-- ^Dense matrix-vector product. Arguments computed in parallel inside the IO monad.
appPIO m v = hv2repa <$> (H.app <$> repa2hmPIO m <*> repa2hvPIO v)

-- Dense matrix-matrix product

mul :: Numeric t => Array F DIM2 t -> Array F DIM2 t -> Array F DIM2 t
-- ^Dense matrix-matrix product.
mul m n = hm2repa $ repa2hm m `H.mul` repa2hm n

mulS :: Numeric t => Array D DIM2 t -> Array D DIM2 t -> Array F DIM2 t
-- ^Dense matrix-matrix product. Arguments computed sequentially.
mulS m n = hm2repa $ repa2hmS m `H.mul` repa2hmS n

mulSIO :: Numeric t => Array D DIM2 t -> Array D DIM2 t -> IO (Array F DIM2 t)
-- ^Dense matrix-matrix product. Arguments computed sequentially inside the IO monad
mulSIO m n = hm2repa <$> (H.mul <$> repa2hmSIO m <*> repa2hmSIO n)

mulP :: (Numeric t, Monad m) => Array D DIM2 t -> Array D DIM2 t -> m (Array F DIM2 t)
-- ^Dense matrix-matrix product. Arguments computed in parallel.
mulP m n = hm2repa <$> (H.mul <$> repa2hmP m <*> repa2hmP n)

mulPIO :: Numeric t => Array D DIM2 t -> Array D DIM2 t -> IO (Array F DIM2 t)
-- ^Dense matrix-matrix product. Arguments computed in parallel inside the IO monad
mulPIO m n = hm2repa <$> (H.mul <$> repa2hmPIO m <*> repa2hmPIO n)

-- Outer product of two vectors

outer :: (Product t, Numeric t) => Array F DIM1 t -> Array F DIM1 t -> Array F DIM2 t
-- |Outer product of two vectors.
outer v u = hm2repa $ repa2hv v `H.outer` repa2hv u

outerS :: (Product t, Numeric t) => Array D DIM1 t -> Array D DIM1 t -> Array F DIM2 t
-- |Outer product of two vectors. Arguments computed sequentially.
outerS v u = hm2repa $ repa2hvS v `H.outer` repa2hvS u

outerSIO :: (Product t, Numeric t) => Array D DIM1 t -> Array D DIM1 t -> IO (Array F DIM2 t)
-- |Outer product of two vectors. Arguments computed sequentially inside the IO monad.
outerSIO v u = hm2repa <$> (H.outer <$> repa2hvSIO v <*> repa2hvSIO u)

outerP :: (Product t, Numeric t, Monad m) => Array D DIM1 t -> Array D DIM1 t -> m (Array F DIM2 t)
-- |Outer product of two vectors. Arguments computed in parallel.
outerP v u = hm2repa <$> (H.outer <$> repa2hvP v <*> repa2hvP u)
outerPIO :: (Product t, Numeric t) => Array D DIM1 t -> Array D DIM1 t -> IO (Array F DIM2 t)
-- |Outer product of two vectors. Arguments computed in parallel inside the IO monad.
outerPIO v u = hm2repa <$> (H.outer <$> repa2hvPIO v <*> repa2hvPIO u)

-- Kronecker product of two matrices

kronecker :: (Product t, Numeric t) => Array F DIM2 t -> Array F DIM2 t -> Array F DIM2 t
-- ^Kronecker product of two matrices.
kronecker m n = hm2repa $ repa2hm m `H.kronecker` repa2hm n

kroneckerS :: (Product t, Numeric t) => Array D DIM2 t -> Array D DIM2 t -> Array F DIM2 t
-- ^Kronecker product of two matrices. Arguments computed sequentially.
kroneckerS m n = hm2repa $ repa2hmS m `H.kronecker` repa2hmS n

kroneckerSIO :: (Product t, Numeric t) => Array D DIM2 t -> Array D DIM2 t -> IO (Array F DIM2 t)
-- ^Kronecker product of two matrices. Arguments computed sequentially inside the IO monad.
kroneckerSIO m n = hm2repa <$> (H.kronecker <$> repa2hmSIO m <*> repa2hmSIO n)

kroneckerP :: (Product t, Numeric t, Monad m) => Array D DIM2 t -> Array D DIM2 t -> m (Array F DIM2 t)
-- ^Kronecker product of two matrices. Arguments computed in parallel.
kroneckerP m n = hm2repa <$> (H.kronecker <$> repa2hmP m <*> repa2hmP n)

kroneckerPIO :: (Product t, Numeric t) => Array D DIM2 t -> Array D DIM2 t -> IO (Array F DIM2 t)
-- ^Kronecker product of two matrices. Arguments computed in parallel inside the IO monad.
kroneckerPIO m n = hm2repa <$> (H.kronecker <$> repa2hmPIO m <*> repa2hmPIO n)

-- Cross product

cross :: Array F DIM1 Double -> Array F DIM1 Double -> Array F DIM1 Double
-- ^Vector cross product.
cross v u = hv2repa $ repa2hv v `H.cross` repa2hv u

crossS :: Array D DIM1 Double -> Array D DIM1 Double -> Array F DIM1 Double
-- ^Vector cross product. Arguments computed sequentially.
crossS v u = hv2repa $ repa2hvS v `H.cross` repa2hvS u

crossSIO :: Array D DIM1 Double -> Array D DIM1 Double -> IO (Array F DIM1 Double)
-- ^Vector cross product. Arguments computed sequentially inside the IO monad.
crossSIO v u = hv2repa <$> (H.cross <$> repa2hvSIO v <*> repa2hvSIO u)

-- Sum of elements

sumElements :: (Numeric t, HShape sh, Container (HType sh) t) => Array F sh t -> t
-- ^Sum elements of a matrix or a vector.
sumElements = H.sumElements . fromRepa

sumElementsS :: (Numeric t, HShape sh, Container (HType sh) t) => Array D sh t -> t
-- ^Sum elements of a matrix or a vector. Argument computed sequentially.
sumElementsS = H.sumElements . fromRepaS

sumElementsSIO :: (Numeric t, HShape sh, Container (HType sh) t) => Array D sh t -> IO t
-- ^Sum elements of a matrix or a vector. Argument computed sequentially in the IO monad.
sumElementsSIO = fmap H.sumElements . fromRepaSIO

sumElementsP :: (Numeric t, HShape sh, Container (HType sh) t, Monad m) => Array D sh t -> m t
-- ^Sum elements of a matrix or a vector. Argument computed in parallel.
sumElementsP = fmap H.sumElements . fromRepaP

sumElementsPIO :: (Numeric t, HShape sh, Container (HType sh) t) => Array D sh t -> IO t
-- ^Sum elements of a matrix or a vector. Argument computed in parallel in the IO monad.
sumElementsPIO = fmap H.sumElements . fromRepaPIO

-- Product of elements

prodElements :: (Numeric t, HShape sh, Container (HType sh) t) => Array F sh t -> t
-- ^Multiply elements of a matrix or a vector.
prodElements = H.prodElements . fromRepa

prodElementsS :: (Numeric t, HShape sh, Container (HType sh) t) => Array D sh t -> t
-- ^Multiply elements of a matrix or a vector. Argument computed sequentially.
prodElementsS = H.prodElements . fromRepaS

prodElementsSIO :: (Numeric t, HShape sh, Container (HType sh) t) => Array D sh t -> IO t
-- ^Multiply elements of a matrix or a vector. Argument computed sequentially inside the IO monad.
prodElementsSIO = fmap H.prodElements . fromRepaSIO

prodElementsP :: (Numeric t, HShape sh, Container (HType sh) t, Monad m) => Array D sh t -> m t
-- ^Multiply elements of a matrix or a vector. Argument computed in parallel.
prodElementsP = fmap H.prodElements . fromRepaP

prodElementsPIO :: (Numeric t, HShape sh, Container (HType sh) t) => Array D sh t -> IO t
-- ^Multiply elements of a matrix or a vector. Argument computed in parallel inside the IO monad.
prodElementsPIO = fmap H.prodElements . fromRepaPIO

-- Linear systems.

(<\>) :: (Field t, Numeric t, HShape sh, LSDiv (HType sh)) => Array F DIM2 t -> Array F sh t -> Array F sh t
-- ^Infix alias for 'solve'.
(<\>) = solve

solve :: (Field t, Numeric t, HShape sh, LSDiv (HType sh)) => Array F DIM2 t -> Array F sh t -> Array F sh t
-- ^Least squares solution of a linear system, similar to the \ operator of Matlab/Octave (based on linearSolveSD).
solve m n = toRepa $ repa2hm m H.<\> fromRepa n

solveS :: (Field t, Numeric t, HShape sh, LSDiv (HType sh)) => Array D DIM2 t -> Array D sh t -> Array F sh t
-- ^Least squares solution of a linear system, similar to the \ operator of Matlab/Octave (based on linearSolveSD). Arguments are computed sequentially.
solveS m n = toRepa $ repa2hmS m H.<\> fromRepaS n

solveSIO :: (Field t, Numeric t, HShape sh, LSDiv (HType sh)) => Array D DIM2 t -> Array D sh t -> IO (Array F sh t)
-- ^Least squares solution of a linear system, similar to the \ operator of Matlab/Octave (based on linearSolveSD). Arguments are computed sequentially inside the IO monad.
solveSIO m n = toRepa <$> ((H.<\>) <$> repa2hmSIO m <*> fromRepaSIO n)

solveP :: (Field t, Numeric t, HShape sh, LSDiv (HType sh), Monad m) => Array D DIM2 t -> Array D sh t -> m (Array F sh t)
-- ^Least squares solution of a linear system, similar to the \ operator of Matlab/Octave (based on linearSolveSD). Arguments are computed in parallel.
solveP m n = toRepa <$> ((H.<\>) <$> repa2hmP m <*> fromRepaP n)

solvePIO :: (Field t, Numeric t, HShape sh, LSDiv (HType sh)) => Array D DIM2 t -> Array D sh t -> IO (Array F sh t)
-- ^Least squares solution of a linear system, similar to the \ operator of Matlab/Octave (based on linearSolveSD). Arguments are computed in parallel inside the IO monad.
solvePIO m n = toRepa <$> ((H.<\>) <$> repa2hmPIO m <*> fromRepaPIO n)


linearSolve :: (Field t, Numeric t) => Array F DIM2 t -> Array F DIM2 t -> Maybe (Array F DIM2 t)
-- ^Solve a linear system (for square coefficient matrix and several right hand sides) using the LU decomposition, returning Nothing for a singular system. For underconstrained or overconstrained systems use 'linearSolveLS' or 'linearSolveSVD'.
linearSolve m n = hm2repa <$> H.linearSolve (repa2hm m) (repa2hm n)

linearSolveS :: (Field t, Numeric t) => Array D DIM2 t -> Array D DIM2 t -> Maybe (Array F DIM2 t)
-- ^Solve a linear system using the LU decomposition. Arguments computed sequentially.
linearSolveS m n = hm2repa <$> H.linearSolve (repa2hmS m) (repa2hmS n)

linearSolveP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> Array D DIM2 t -> m (Maybe (Array F DIM2 t))
-- ^Solve a linear system using the LU decomposition. Arguments computed in parallel.
linearSolveP m n = (hm2repa <$>) <$> (H.linearSolve <$> repa2hmP m <*> repa2hmP n)

linearSolveSIO :: (Field t, Numeric t) => Array D DIM2 t -> Array D DIM2 t -> IO (Maybe (Array F DIM2 t))
-- ^Solve a linear system using the LU decomposition. Arguments computed sequentially inside the IO monad.
linearSolveSIO m n = (hm2repa <$>) <$> (H.linearSolve <$> repa2hmSIO m <*> repa2hmSIO n)

linearSolvePIO :: (Field t, Numeric t) => Array D DIM2 t -> Array D DIM2 t -> IO (Maybe (Array F DIM2 t))
-- ^Solve a linear system using the LU decomposition. Arguments computed in parallel inside the IO monad.
linearSolvePIO m n = (hm2repa <$>) <$> (H.linearSolve <$> repa2hmPIO m <*> repa2hmPIO n)


linearSolveLS :: (Field t, Numeric t) => Array F DIM2 t -> Array F DIM2 t -> Array F DIM2 t
-- ^Least squared error solution of an overcompensated system, or the minimum norm solution of an undercompensated system. For rank-deficient systems use 'linearSolveSVD'.
linearSolveLS m n = hm2repa $ H.linearSolveLS (repa2hm m) (repa2hm n)

linearSolveLS_S :: (Field t, Numeric t) => Array D DIM2 t -> Array D DIM2 t -> Array F DIM2 t 
linearSolveLS_S m n = hm2repa $ H.linearSolveLS (repa2hmS m) (repa2hmS n)

linearSolveLS_SIO :: (Field t, Numeric t) => Array D DIM2 t -> Array D DIM2 t -> IO (Array F DIM2 t)
linearSolveLS_SIO m n = hm2repa <$> (H.linearSolveLS <$> repa2hmSIO m <*> repa2hmSIO n)

linearSolveLS_P :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> Array D DIM2 t -> m (Array F DIM2 t)
linearSolveLS_P m n = hm2repa <$> (H.linearSolveLS <$> repa2hmP m <*> repa2hmP n)

linearSolveLS_PIO :: (Field t, Numeric t) => Array D DIM2 t -> Array D DIM2 t -> IO (Array F DIM2 t)
linearSolveLS_PIO m n = hm2repa <$> (H.linearSolveLS <$> repa2hmPIO m <*> repa2hmPIO n)


linearSolveSVD :: (Field t, Numeric t) => Array F DIM2 t -> Array F DIM2 t -> Array F DIM2 t
-- ^Minimum norm solution of a general linear least squares problem Ax=b using the SVD. Admits rank-deficient systems but is slower than 'linearSolveLS'. The effective rank of A is determined by treating as zero those singular values which are less than eps times the largest singular value.
linearSolveSVD m n = hm2repa $ H.linearSolveSVD (repa2hm m) (repa2hm n)

linearSolveSVD_S :: (Field t, Numeric t) => Array D DIM2 t -> Array D DIM2 t -> Array F DIM2 t 
linearSolveSVD_S m n = hm2repa $ H.linearSolveSVD (repa2hmS m) (repa2hmS n)

linearSolveSVD_SIO :: (Field t, Numeric t) => Array D DIM2 t -> Array D DIM2 t -> IO (Array F DIM2 t)
linearSolveSVD_SIO m n = hm2repa <$> (H.linearSolveSVD <$> repa2hmSIO m <*> repa2hmSIO n)

linearSolveSVD_P :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> Array D DIM2 t -> m (Array F DIM2 t)
linearSolveSVD_P m n = hm2repa <$> (H.linearSolveSVD <$> repa2hmP m <*> repa2hmP n)

linearSolveSVD_PIO :: (Field t, Numeric t) => Array D DIM2 t -> Array D DIM2 t -> IO (Array F DIM2 t)
linearSolveSVD_PIO m n = hm2repa <$> (H.linearSolveLS <$> repa2hmPIO m <*> repa2hmPIO n)


luSolve :: (Field t, Numeric t) => (Array F DIM2 t, [Int]) -> Array F DIM2 t -> Array F DIM2 t
-- ^Solution of a linear system (for several right hand sides) from the precomputed LU factorization obtained by 'luPacked'.
luSolve (lu, l) m = hm2repa $ H.luSolve (repa2hm lu, l) (repa2hm m)

luSolveS :: (Field t, Numeric t) => (Array F DIM2 t, [Int]) -> Array D DIM2 t -> Array F DIM2 t
luSolveS (lu, l) m = hm2repa $ H.luSolve (repa2hm lu, l) (repa2hmS m)

luSolveSIO :: (Field t, Numeric t) => (Array F DIM2 t, [Int]) -> Array D DIM2 t -> IO (Array F DIM2 t)
luSolveSIO (lu, l) m = hm2repa . H.luSolve (repa2hm lu, l) <$> repa2hmSIO m

luSolveP :: (Field t, Numeric t, Monad m) => (Array F DIM2 t, [Int]) -> Array D DIM2 t -> m (Array F DIM2 t)
luSolveP (lu, l) m = hm2repa . H.luSolve (repa2hm lu, l) <$> repa2hmP m

luSolvePIO :: (Field t, Numeric t) => (Array F DIM2 t, [Int]) -> Array D DIM2 t -> IO (Array F DIM2 t)
luSolvePIO (lu, l) m = hm2repa . H.luSolve (repa2hm lu, l) <$> repa2hmPIO m


cholSolve :: (Field t, Numeric t) => Array F DIM2 t -> Array F DIM2 t -> Array F DIM2 t
-- ^Solve a symmetric or Herimitian positive definite linear system using a precomputed Cholesky decomposition obtained by 'chol'.
cholSolve ch m = hm2repa $ H.cholSolve (repa2hm ch) (repa2hm m)

cholSolveS :: (Field t, Numeric t) => Array F DIM2 t -> Array D DIM2 t -> Array F DIM2 t
cholSolveS ch m = hm2repa $ H.cholSolve (repa2hm ch) (repa2hmS m)

cholSolveSIO :: (Field t, Numeric t) => Array F DIM2 t -> Array D DIM2 t -> IO (Array F DIM2 t)
cholSolveSIO ch m = hm2repa . H.cholSolve (repa2hm ch) <$> repa2hmSIO m

cholSolveP :: (Field t, Numeric t, Monad m) => Array F DIM2 t -> Array D DIM2 t -> m (Array F DIM2 t)
cholSolveP ch m = hm2repa . H.cholSolve (repa2hm ch) <$> repa2hmP m

cholSolvePIO :: (Field t, Numeric t) => Array F DIM2 t -> Array D DIM2 t -> IO (Array F DIM2 t)
cholSolvePIO ch m = hm2repa . H.cholSolve (repa2hm ch) <$> repa2hmPIO m

-- Inverse and pseudoinverse

inv :: (Field t, Numeric t) => Array F DIM2 t -> Array F DIM2 t
-- ^Inverse of a square matrix.
inv = hm2repa . H.inv . repa2hm

invS :: (Field t, Numeric t) => Array D DIM2 t -> Array F DIM2 t
invS = hm2repa . H.inv . repa2hmS

invSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t)
invSIO = fmap (hm2repa . H.inv) . repa2hmSIO

invP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM2 t)
invP = fmap (hm2repa . H.inv) . repa2hmP

invPIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t)
invPIO = fmap (hm2repa . H.inv) . repa2hmPIO

pinv :: (Field t, Numeric t) => Array F DIM2 t -> Array F DIM2 t
-- ^Pseudoinverse of a general matrix, with default tolerance ('pinvTol' 1, similar to GNU-Octave)
pinv = hm2repa . H.pinv . repa2hm

pinvS :: (Field t, Numeric t) => Array D DIM2 t -> Array F DIM2 t
pinvS = hm2repa . H.pinv . repa2hmS

pinvSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t)
pinvSIO = fmap (hm2repa . H.pinv) . repa2hmSIO

pinvP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM2 t)
pinvP = fmap (hm2repa . H.pinv) . repa2hmP

pinvPIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t)
pinvPIO = fmap (hm2repa . H.pinv) . repa2hmPIO

pinvTol :: (Field t, Numeric t) => Double -> Array F DIM2 t -> Array F DIM2 t
-- ^pinvTol r computes the pseudoinverse of a matrix with tolerance tol=r*g*eps*(max rows cols), where g is the greatest singular value.
pinvTol r = hm2repa . H.pinvTol r . repa2hm

pinvTolS :: (Field t, Numeric t) => Double -> Array D DIM2 t -> Array F DIM2 t
pinvTolS r = hm2repa . H.pinvTol r . repa2hmS

pinvTolSIO :: (Field t, Numeric t) => Double -> Array D DIM2 t -> IO (Array F DIM2 t)
pinvTolSIO r= fmap (hm2repa . H.pinvTol r) . repa2hmSIO

pinvTolP :: (Field t, Numeric t, Monad m) => Double -> Array D DIM2 t -> m (Array F DIM2 t)
pinvTolP r = fmap (hm2repa . H.pinvTol r) . repa2hmP

pinvTolPIO :: (Field t, Numeric t) => Double -> Array D DIM2 t -> IO (Array F DIM2 t)
pinvTolPIO r = fmap (hm2repa . H.pinvTol r) . repa2hmPIO

-- Determinant and rank

rcond :: (Field t, Numeric t) => Array F DIM2 t -> Double
-- ^Reciprocal of the 2-norm condition number of a matrix, computed from the singular values.
rcond = H.rcond . repa2hm

rcondS :: (Field t, Numeric t) => Array D DIM2 t -> Double
rcondS = H.rcond . repa2hmS

rcondSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO Double
rcondSIO = fmap H.rcond . repa2hmSIO

rcondP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m Double
rcondP = fmap H.rcond . repa2hmP

rcondPIO :: (Field t, Numeric t) => Array D DIM2 t -> IO Double
rcondPIO = fmap H.rcond . repa2hmPIO

rank :: (Field t, Numeric t) => Array F DIM2 t -> Int
-- ^Number of linearly independent rows or columns. See also 'ranksv'.
rank = H.rank . repa2hm

rankS :: (Field t, Numeric t) => Array D DIM2 t -> Int
rankS = H.rank . repa2hmS

rankSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO Int
rankSIO = fmap H.rank . repa2hmSIO

rankP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m Int
rankP = fmap H.rank . repa2hmP

rankPIO :: (Field t, Numeric t) => Array D DIM2 t -> IO Int
rankPIO = fmap H.rank . repa2hmPIO

det :: (Field t, Numeric t) => Array F DIM2 t -> t
-- ^Determinant of a square matrix. To avoid possible overflow or underflow use 'invlndet'.
det = H.det . repa2hm

detS :: (Field t, Numeric t) => Array D DIM2 t -> t
detS = H.det . repa2hmS

detSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO t
detSIO = fmap H.det . repa2hmSIO

detP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m t
detP = fmap H.det . repa2hmP

detPIO :: (Field t, Numeric t) => Array D DIM2 t -> IO t
detPIO = fmap H.det . repa2hmPIO

invlndet :: (Field t, Numeric t) => Array F DIM2 t -> (Array F DIM2 t, (t, t)) -- ^(inverse, (log abs det, sign or phase of det))
-- ^Joint computation of inverse and logarithm of determinant of a square matrix.
invlndet m = let (h, r) = H.invlndet $ repa2hm m in (hm2repa h, r)
  
invlndetS :: (Field t, Numeric t) => Array D DIM2 t -> (Array F DIM2 t, (t, t))
invlndetS m = let (h, r) = H.invlndet $ repa2hmS m in (hm2repa h, r)

invlndetSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, (t, t))
invlndetSIO m = do 
  (h, r) <- H.invlndet <$> repa2hmSIO m 
  return (hm2repa h, r)

invlndetP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM2 t, (t, t))
invlndetP m = do 
  (h, r) <- H.invlndet <$> repa2hmP m 
  return (hm2repa h, r)

invlndetPIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, (t, t))
invlndetPIO m = do 
  (h, r) <- H.invlndet <$> repa2hmPIO m 
  return (hm2repa h, r)

-- Norms

norm_Frob :: (Normed (Vector t), Element t) => Array F DIM2 t -> Double
norm_Frob = H.norm_Frob . repa2hm

norm_FrobS :: (Normed (Vector t), Element t) => Array D DIM2 t -> Double
norm_FrobS = H.norm_Frob . repa2hmS

norm_FrobSIO :: (Normed (Vector t), Element t) => Array D DIM2 t -> IO Double
norm_FrobSIO = fmap H.norm_Frob . repa2hmSIO

norm_FrobP :: (Normed (Vector t), Element t, Monad m) => Array D DIM2 t -> m Double
norm_FrobP = fmap H.norm_Frob . repa2hmP

norm_FrobPIO :: (Normed (Vector t), Element t) => Array D DIM2 t -> IO Double
norm_FrobPIO = fmap H.norm_Frob . repa2hmPIO


norm_nuclear :: (Field t, Numeric t) => Array F DIM2 t -> Double
norm_nuclear = H.norm_nuclear . repa2hm

norm_nuclearS :: (Field t, Numeric t) => Array D DIM2 t -> Double
norm_nuclearS = H.norm_nuclear . repa2hmS

norm_nuclearSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO Double
norm_nuclearSIO = fmap H.norm_nuclear . repa2hmSIO

norm_nuclearP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m Double
norm_nuclearP = fmap H.norm_nuclear . repa2hmP

norm_nuclearPIO :: (Field t, Numeric t) => Array D DIM2 t -> IO Double
norm_nuclearPIO = fmap H.norm_nuclear . repa2hmPIO

-- Nullspace and range

orth :: (Field t, Numeric t) => Array F DIM2 t -> Array F DIM2 t
-- ^An orthonormal basis of the range space of a matrix. See also 'orthSVD'.
orth = hm2repa . H.orth . repa2hm

orthS :: (Field t, Numeric t) => Array D DIM2 t -> Array F DIM2 t
orthS = hm2repa . H.orth . repa2hmS

orthSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t)
orthSIO = fmap (hm2repa . H.orth) . repa2hmSIO

orthP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM2 t)
orthP = fmap (hm2repa . H.orth) . repa2hmP

orthPIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t)
orthPIO = fmap (hm2repa . H.orth) . repa2hmPIO

nullspace :: (Field t, Numeric t) => Array F DIM2 t -> Array F DIM2 t
-- ^An orthonormal basis of the null space of a matrix. See also 'nullspaceSVD'.
nullspace = hm2repa . H.nullspace . repa2hm

nullspaceS :: (Field t, Numeric t) => Array D DIM2 t -> Array F DIM2 t
nullspaceS = hm2repa . H.nullspace . repa2hmS

nullspaceSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t)
nullspaceSIO = fmap (hm2repa . H.nullspace) . repa2hmSIO

nullspaceP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM2 t)
nullspaceP = fmap (hm2repa . H.nullspace) . repa2hmP

nullspacePIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t)
nullspacePIO = fmap (hm2repa . H.nullspace) . repa2hmPIO

null1 :: Array F DIM2 Double -> Array F DIM1 Double
-- ^Solution of an overconstrained homogenous linear system.
null1 = hv2repa . H.null1 . repa2hm

null1S :: Array D DIM2 Double -> Array F DIM1 Double
null1S = hv2repa . H.null1 . repa2hmS

null1SIO :: Array D DIM2 Double -> IO (Array F DIM1 Double)
null1SIO = fmap (hv2repa . H.null1) . repa2hmSIO

null1P :: Monad m => Array D DIM2 Double -> m (Array F DIM1 Double)
null1P = fmap (hv2repa . H.null1) . repa2hmP

null1PIO :: Array D DIM2 Double -> IO (Array F DIM1 Double)
null1PIO = fmap (hv2repa . H.null1) . repa2hmPIO

null1sym :: Array F DIM2 Double -> Array F DIM1 Double
-- ^Solution of an overconstrained homogenous symmetric linear system.
null1sym = hv2repa . H.null1sym . repa2hm

null1symS :: Array D DIM2 Double -> Array F DIM1 Double
null1symS = hv2repa . H.null1sym . repa2hmS

null1symSIO :: Array D DIM2 Double -> IO (Array F DIM1 Double)
null1symSIO = fmap (hv2repa . H.null1sym) . repa2hmSIO

null1symP :: Monad m => Array D DIM2 Double -> m (Array F DIM1 Double)
null1symP = fmap (hv2repa . H.null1sym) . repa2hmP

null1symPIO :: Array D DIM2 Double -> IO (Array F DIM1 Double)
null1symPIO = fmap (hv2repa . H.null1sym) . repa2hmPIO

-- SVD

svd :: (Field t, Numeric t) => Array F DIM2 t -> (Array F DIM2 t, Array F DIM1 Double, Array F DIM2 t)
-- ^Full singular value decomposition.
svd m = let (u,s,v) = H.svd $ repa2hm m in (hm2repa u, hv2repa s, hm2repa v)

svdS :: (Field t, Numeric t) => Array D DIM2 t -> (Array F DIM2 t, Array F DIM1 Double, Array F DIM2 t)
svdS m = let (u,s,v) = H.svd $ repa2hmS m in (hm2repa u, hv2repa s, hm2repa v)

svdSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, Array F DIM1 Double, Array F DIM2 t)
svdSIO m = do
  (u,s,v) <- H.svd <$> repa2hmSIO m 
  return (hm2repa u, hv2repa s, hm2repa v)

svdP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM2 t, Array F DIM1 Double, Array F DIM2 t)
svdP m = do
  (u,s,v) <- H.svd <$> repa2hmP m 
  return (hm2repa u, hv2repa s, hm2repa v)

svdPIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, Array F DIM1 Double, Array F DIM2 t)
svdPIO m = do
  (u,s,v) <- H.svd <$> repa2hmPIO m 
  return (hm2repa u, hv2repa s, hm2repa v)

thinSVD :: (Field t, Numeric t) => Array F DIM2 t -> (Array F DIM2 t, Array F DIM1 Double, Array F DIM2 t)
-- ^A version of 'svd' which returns only the min (rows m) (cols m) singular vectors of m. (u,s,v) = thinSVD m ==> m == u * diag s * tr v
thinSVD m = let (u,s,v) = H.thinSVD $ repa2hm m in (hm2repa u, hv2repa s, hm2repa v)

thinSVD_S :: (Field t, Numeric t) => Array D DIM2 t -> (Array F DIM2 t, Array F DIM1 Double, Array F DIM2 t)
thinSVD_S m = let (u,s,v) = H.thinSVD $ repa2hmS m in (hm2repa u, hv2repa s, hm2repa v)

thinSVD_SIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, Array F DIM1 Double, Array F DIM2 t)
thinSVD_SIO m = do
  (u,s,v) <- H.thinSVD <$> repa2hmSIO m 
  return (hm2repa u, hv2repa s, hm2repa v)

thinSVD_P :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM2 t, Array F DIM1 Double, Array F DIM2 t)
thinSVD_P m = do
  (u,s,v) <- H.thinSVD <$> repa2hmP m 
  return (hm2repa u, hv2repa s, hm2repa v)

thinSVD_PIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, Array F DIM1 Double, Array F DIM2 t)
thinSVD_PIO m = do
  (u,s,v) <- H.thinSVD <$> repa2hmPIO m 
  return (hm2repa u, hv2repa s, hm2repa v)

compactSVD :: (Field t, Numeric t) => Array F DIM2 t -> (Array F DIM2 t, Array F DIM1 Double, Array F DIM2 t)
-- ^Similar to 'thinSVD', returning only the nonzero singular values and the corresponding singular vectors.
compactSVD m = let (u,s,v) = H.compactSVD $ repa2hm m in (hm2repa u, hv2repa s, hm2repa v)

compactSVD_S :: (Field t, Numeric t) => Array D DIM2 t -> (Array F DIM2 t, Array F DIM1 Double, Array F DIM2 t)
compactSVD_S m = let (u,s,v) = H.compactSVD $ repa2hmS m in (hm2repa u, hv2repa s, hm2repa v)

compactSVD_SIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, Array F DIM1 Double, Array F DIM2 t)
compactSVD_SIO m = do
  (u,s,v) <- H.compactSVD <$> repa2hmSIO m 
  return (hm2repa u, hv2repa s, hm2repa v)

compactSVD_P :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM2 t, Array F DIM1 Double, Array F DIM2 t)
compactSVD_P m = do
  (u,s,v) <- H.compactSVD <$> repa2hmP m 
  return (hm2repa u, hv2repa s, hm2repa v)

compactSVD_PIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, Array F DIM1 Double, Array F DIM2 t)
compactSVD_PIO m = do
  (u,s,v) <- H.compactSVD <$> repa2hmPIO m 
  return (hm2repa u, hv2repa s, hm2repa v)

singularValues :: (Field t, Numeric t) => Array F DIM2 t -> Array F DIM1 Double
-- ^Singular values only.
singularValues = hv2repa . H.singularValues . repa2hm

singularValuesS :: (Field t, Numeric t) => Array D DIM2 t -> Array F DIM1 Double
singularValuesS = hv2repa . H.singularValues . repa2hmS

singularValuesSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM1 Double)
singularValuesSIO = fmap (hv2repa . H.singularValues) . repa2hmSIO

singularValuesP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM1 Double)
singularValuesP = fmap (hv2repa . H.singularValues) . repa2hmP

singularValuesPIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM1 Double)
singularValuesPIO = fmap (hv2repa . H.singularValues) . repa2hmPIO

leftSV :: (Field t, Numeric t) => Array F DIM2 t -> (Array F DIM2 t, Array F DIM1 Double)
-- ^Singular values and all left singular vectors (as columns).
leftSV m = let (u,s) = H.leftSV $ repa2hm m in (hm2repa u, hv2repa s)

leftSV_S :: (Field t, Numeric t) => Array D DIM2 t -> (Array F DIM2 t, Array F DIM1 Double)
leftSV_S m = let (u,s) = H.leftSV $ repa2hmS m in (hm2repa u, hv2repa s)

leftSV_SIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, Array F DIM1 Double)
leftSV_SIO m = do
  (u,s) <- H.leftSV <$> repa2hmSIO m 
  return (hm2repa u, hv2repa s)

leftSV_P :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM2 t, Array F DIM1 Double)
leftSV_P m = do
  (u,s) <- H.leftSV <$> repa2hmP m 
  return (hm2repa u, hv2repa s)

leftSV_PIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, Array F DIM1 Double)
leftSV_PIO m = do
  (u,s) <- H.leftSV <$> repa2hmPIO m 
  return (hm2repa u, hv2repa s)

rightSV :: (Field t, Numeric t) => Array F DIM2 t -> (Array F DIM1 Double, Array F DIM2 t)
-- ^Singular values and all right singular vectors (as columns).
rightSV m = let (s,v) = H.rightSV $ repa2hm m in (hv2repa s, hm2repa v)

rightSV_S :: (Field t, Numeric t) => Array D DIM2 t -> (Array F DIM1 Double, Array F DIM2 t)
rightSV_S m = let (s,v) = H.rightSV $ repa2hmS m in (hv2repa s, hm2repa v)

rightSV_SIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM1 Double, Array F DIM2 t)
rightSV_SIO m = do
  (s,v) <- H.rightSV <$> repa2hmSIO m 
  return (hv2repa s, hm2repa v)

rightSV_P :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM1 Double, Array F DIM2 t)
rightSV_P m = do
  (s,v) <- H.rightSV <$> repa2hmP m 
  return (hv2repa s, hm2repa v)

rightSV_PIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM1 Double, Array F DIM2 t)
rightSV_PIO m = do
  (s,v) <- H.rightSV <$> repa2hmPIO m 
  return (hv2repa s, hm2repa v)

-- Eigensystems

eig :: (Field t, Numeric t) => Array F DIM2 t -> (Array F DIM1 (Complex Double), Array F DIM2 (Complex Double))
-- ^Eigenvalues (not ordered) and eigenvectors (as columns) of a general square matrix. (s,v) = eig m ==> m * v = v == v <> diag s
eig m = let (s,v) = H.eig $ repa2hm m in (hv2repa s, hm2repa v)

eigS :: (Field t, Numeric t) => Array D DIM2 t -> (Array F DIM1 (Complex Double), Array F DIM2 (Complex Double))
eigS m = let (s,v) = H.eig $ repa2hmS m in (hv2repa s, hm2repa v)

eigSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM1 (Complex Double), Array F DIM2 (Complex Double))
eigSIO m = do
  (s,v) <- H.eig <$> repa2hmSIO m 
  return (hv2repa s, hm2repa v)

eigP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM1 (Complex Double), Array F DIM2 (Complex Double))
eigP m = do
  (s,v) <- H.eig <$> repa2hmP m 
  return (hv2repa s, hm2repa v)

eigPIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM1 (Complex Double), Array F DIM2 (Complex Double))
eigPIO m = do
  (s,v) <- H.eig <$> repa2hmPIO m 
  return (hv2repa s, hm2repa v)

eigSH :: (Field t, Numeric t) => Array F DIM2 t -> (Array F DIM1 Double, Array F DIM2 t)
-- ^Eigenvalues and eigenvectors (as columns) of a complex hermitian or a real symmetric matrix, in descending order.
eigSH m = let (s,v) = H.eigSH $ repa2hm m in (hv2repa s, hm2repa v)

eigSH_S :: (Field t, Numeric t) => Array D DIM2 t -> (Array F DIM1 Double, Array F DIM2 t)
eigSH_S m = let (s,v) = H.eigSH $ repa2hmS m in (hv2repa s, hm2repa v)

eigSH_SIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM1 Double, Array F DIM2 t)
eigSH_SIO m = do
  (s,v) <- H.eigSH <$> repa2hmSIO m 
  return (hv2repa s, hm2repa v)

eigSH_P :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM1 Double, Array F DIM2 t)
eigSH_P m = do
  (s,v) <- H.eigSH <$> repa2hmP m 
  return (hv2repa s, hm2repa v)

eigSH_PIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM1 Double, Array F DIM2 t)
eigSH_PIO m = do
  (s,v) <- H.eigSH <$> repa2hmPIO m 
  return (hv2repa s, hm2repa v)

eigSH' :: (Field t, Numeric t) => Array F DIM2 t -> (Array F DIM1 Double, Array F DIM2 t)
-- ^Similar to 'eigSH' without checking that the input matrix is hermitian or symmetric. It works with the upper triangular part.
eigSH' m = let (s,v) = H.eigSH' $ repa2hm m in (hv2repa s, hm2repa v)

eigSH'S :: (Field t, Numeric t) => Array D DIM2 t -> (Array F DIM1 Double, Array F DIM2 t)
eigSH'S m = let (s,v) = H.eigSH' $ repa2hmS m in (hv2repa s, hm2repa v)

eigSH'SIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM1 Double, Array F DIM2 t)
eigSH'SIO m = do
  (s,v) <- H.eigSH' <$> repa2hmSIO m 
  return (hv2repa s, hm2repa v)

eigSH'P :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM1 Double, Array F DIM2 t)
eigSH'P m = do
  (s,v) <- H.eigSH' <$> repa2hmP m 
  return (hv2repa s, hm2repa v)

eigSH'PIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM1 Double, Array F DIM2 t)
eigSH'PIO m = do
  (s,v) <- H.eigSH' <$> repa2hmPIO m 
  return (hv2repa s, hm2repa v)

eigenvalues :: (Field t, Numeric t) => Array F DIM2 t -> Array F DIM1 (Complex Double)
-- ^Eigenvalues (not ordered) of a general square matrix.
eigenvalues = hv2repa . H.eigenvalues . repa2hm

eigenvaluesS :: (Field t, Numeric t) => Array D DIM2 t -> Array F DIM1 (Complex Double)
eigenvaluesS = hv2repa . H.eigenvalues . repa2hmS

eigenvaluesSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM1 (Complex Double))
eigenvaluesSIO = fmap (hv2repa . H.eigenvalues) . repa2hmSIO

eigenvaluesP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM1 (Complex Double))
eigenvaluesP = fmap (hv2repa . H.eigenvalues) . repa2hmP

eigenvaluesPIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM1 (Complex Double))
eigenvaluesPIO = fmap (hv2repa . H.eigenvalues) . repa2hmPIO

eigenvaluesSH :: (Field t, Numeric t) => Array F DIM2 t -> Array F DIM1 Double
-- ^Eigenvalues (in descending order) of a complex hermitian or real symmetric matrix.
eigenvaluesSH = hv2repa . H.eigenvaluesSH . repa2hm

eigenvaluesSH_S :: (Field t, Numeric t) => Array D DIM2 t -> Array F DIM1 Double
eigenvaluesSH_S = hv2repa . H.eigenvaluesSH . repa2hmS

eigenvaluesSH_SIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM1 Double)
eigenvaluesSH_SIO = fmap (hv2repa . H.eigenvaluesSH) . repa2hmSIO

eigenvaluesSH_P :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM1 Double)
eigenvaluesSH_P = fmap (hv2repa . H.eigenvaluesSH) . repa2hmP

eigenvaluesSH_PIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM1 Double)
eigenvaluesSH_PIO = fmap (hv2repa . H.eigenvaluesSH) . repa2hmPIO

eigenvaluesSH' :: (Field t, Numeric t) => Array F DIM2 t -> Array F DIM1 Double
-- ^Similar to 'eigenvaluesSH' without checking that the input matrix is hermitian or symmetric. It works with the upper triangular part.
eigenvaluesSH' = hv2repa . H.eigenvaluesSH' . repa2hm

eigenvaluesSH'S :: (Field t, Numeric t) => Array D DIM2 t -> Array F DIM1 Double
eigenvaluesSH'S = hv2repa . H.eigenvaluesSH' . repa2hmS

eigenvaluesSH'SIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM1 Double)
eigenvaluesSH'SIO = fmap (hv2repa . H.eigenvaluesSH') . repa2hmSIO

eigenvaluesSH'P :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM1 Double)
eigenvaluesSH'P = fmap (hv2repa . H.eigenvaluesSH') . repa2hmP

eigenvaluesSH'PIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM1 Double)
eigenvaluesSH'PIO = fmap (hv2repa . H.eigenvaluesSH') . repa2hmPIO

geigSH' :: (Field t, Numeric t) => Array F DIM2 t -> Array F DIM2 t -> (Array F DIM1 Double, Array F DIM2 t)
-- ^Generalized symmetric positive definite eigensystem Av = IBv, for A and B symmetric, B positive definite (conditions not checked).
geigSH' a b = let (s,v) = H.geigSH' (repa2hm a) (repa2hm b) in (hv2repa s, hm2repa v)

geigSH'S :: (Field t, Numeric t) => Array D DIM2 t -> Array D DIM2 t -> (Array F DIM1 Double, Array F DIM2 t)
geigSH'S a b = let (s,v) = H.geigSH' (repa2hmS a) (repa2hmS b) in (hv2repa s, hm2repa v)

geigSH'SIO :: (Field t, Numeric t) => Array D DIM2 t -> Array D DIM2 t -> IO (Array F DIM1 Double, Array F DIM2 t)
geigSH'SIO a b = do
  (s,v) <- H.geigSH' <$> repa2hmSIO a <*> repa2hmSIO b
  return (hv2repa s, hm2repa v)

geigSH'P :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> Array D DIM2 t -> m (Array F DIM1 Double, Array F DIM2 t)
geigSH'P a b = do
  (s,v) <- H.geigSH' <$> repa2hmP a <*> repa2hmP b
  return (hv2repa s, hm2repa v)

geigSH'PIO :: (Field t, Numeric t) => Array D DIM2 t -> Array D DIM2 t -> IO (Array F DIM1 Double, Array F DIM2 t)
geigSH'PIO a b = do
  (s,v) <- H.geigSH' <$> repa2hmPIO a <*> repa2hmPIO b
  return (hv2repa s, hm2repa v)

-- QR

qr :: (Field t, Numeric t) => Array F DIM2 t -> (Array F DIM2 t, Array F DIM2 t)
-- ^QR factorization. (q,r) = qr m ==> m = q * r where q is unitary and r is upper triangular.
qr m = let (q,r) = H.qr $ repa2hm m in (hm2repa q, hm2repa r)

qrS :: (Field t, Numeric t) => Array D DIM2 t -> (Array F DIM2 t, Array F DIM2 t)
qrS m = let (q,r) = H.qr $ repa2hmS m in (hm2repa q, hm2repa r)

qrSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, Array F DIM2 t)
qrSIO m = do
  (q,r) <- H.qr <$> repa2hmSIO m 
  return (hm2repa q, hm2repa r)

qrP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM2 t, Array F DIM2 t)
qrP m = do
  (q,r) <- H.qr <$> repa2hmP m 
  return (hm2repa q, hm2repa r)

qrPIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, Array F DIM2 t)
qrPIO m = do
  (q,r) <- H.qr <$> repa2hmPIO m 
  return (hm2repa q, hm2repa r)

rq :: (Field t, Numeric t) => Array F DIM2 t -> (Array F DIM2 t, Array F DIM2 t)
-- ^RQ factorization. (r,q) = rq m ==> m = r * q where q is unitary and r is upper triangular.
rq m = let (r,q) = H.rq $ repa2hm m in (hm2repa r, hm2repa q)

rqS :: (Field t, Numeric t) => Array D DIM2 t -> (Array F DIM2 t, Array F DIM2 t)
rqS m = let (r,q) = H.rq $ repa2hmS m in (hm2repa r, hm2repa q)

rqSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, Array F DIM2 t)
rqSIO m = do 
  (r,q) <- H.rq <$> repa2hmSIO m 
  return (hm2repa r, hm2repa q)

rqP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM2 t, Array F DIM2 t)
rqP m = do 
  (r,q) <- H.rq <$> repa2hmP m 
  return (hm2repa r, hm2repa q)

rqPIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, Array F DIM2 t)
rqPIO m = do 
  (r,q) <- H.rq <$> repa2hmPIO m 
  return (hm2repa r, hm2repa q)

qrRaw :: (Field t, Numeric t) => Array F DIM2 t -> (Array F DIM2 t, Array F DIM1 t)
qrRaw m = let (n,v) = H.qrRaw $ repa2hm m in (hm2repa n, hv2repa v)

qrRawS :: (Field t, Numeric t) => Array D DIM2 t -> (Array F DIM2 t, Array F DIM1 t)
qrRawS m = let (n,v) = H.qrRaw $ repa2hmS m in (hm2repa n, hv2repa v)

qrRawSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, Array F DIM1 t)
qrRawSIO m = do
  (n,v) <- H.qrRaw <$> repa2hmSIO m 
  return (hm2repa n, hv2repa v)

qrRawP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM2 t, Array F DIM1 t)
qrRawP m = do
  (n,v) <- H.qrRaw <$> repa2hmP m 
  return (hm2repa n, hv2repa v)

qrRawPIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, Array F DIM1 t)
qrRawPIO m = do
  (n,v) <- H.qrRaw <$> repa2hmPIO m 
  return (hm2repa n, hv2repa v)

qrgr :: (Field t, Numeric t) => Int -> (Array F DIM2 t, Array F DIM1 t) -> Array F DIM2 t
-- ^Generate a matrix with k othogonal columns from the output of 'qrRaw'.
qrgr k (m,v) = hm2repa $ H.qrgr k (repa2hm m, repa2hv v)

-- Cholesky

chol :: (Field t, Numeric t) => Array F DIM2 t -> Array F DIM2 t
-- ^Cholesky factorization of a positive definite hermitian or symmetric matrix. c = chol m ==> m == c' * c where c is upper triangular.
chol = hm2repa . H.chol . repa2hm

cholS :: (Field t, Numeric t) => Array D DIM2 t -> Array F DIM2 t
cholS = hm2repa . H.chol . repa2hmS

cholSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t)
cholSIO = fmap (hm2repa . H.chol) . repa2hmSIO

cholP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM2 t)
cholP = fmap (hm2repa . H.chol) . repa2hmP

cholPIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t)
cholPIO = fmap (hm2repa . H.chol) . repa2hmPIO

chol' :: (Field t, Numeric t) => Array F DIM2 t -> Array F DIM2 t
-- ^Similar to 'chol' without checking that the input matrix is hermitian or symmetric. It works with the upper triangular part.
chol' = hm2repa . H.cholSH . repa2hm

chol'S :: (Field t, Numeric t) => Array D DIM2 t -> Array F DIM2 t
chol'S = hm2repa . H.cholSH . repa2hmS

chol'SIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t)
chol'SIO = fmap (hm2repa . H.cholSH) . repa2hmSIO

chol'P :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM2 t)
chol'P = fmap (hm2repa . H.cholSH) . repa2hmP

chol'PIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t)
chol'PIO = fmap (hm2repa . H.cholSH) . repa2hmPIO

-- Hessenberg

hess :: (Field t, Numeric t) => Array F DIM2 t -> (Array F DIM2 t, Array F DIM2 t)
-- ^Hessenberg factorization. (p,h) == hess m ==> p * h * p' where p is unitary and h is in upper Hessenberg form (zero entries below the first subdiagonal).
hess m = let (p,h) = H.hess $ repa2hm m in (hm2repa p, hm2repa h)

hessS :: (Field t, Numeric t) => Array D DIM2 t -> (Array F DIM2 t, Array F DIM2 t)
hessS m = let (p,h) = H.hess $ repa2hmS m in (hm2repa p, hm2repa h)

hessSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, Array F DIM2 t)
hessSIO m = do 
  (p,h) <- H.hess <$> repa2hmSIO m 
  return (hm2repa p, hm2repa h)

hessP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM2 t, Array F DIM2 t)
hessP m = do 
  (p,h) <- H.hess <$> repa2hmP m 
  return (hm2repa p, hm2repa h)

hessPIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, Array F DIM2 t)
hessPIO m = do 
  (p,h) <- H.hess <$> repa2hmPIO m 
  return (hm2repa p, hm2repa h)

-- Schur
schur :: (Field t, Numeric t) => Array F DIM2 t -> (Array F DIM2 t, Array F DIM2 t)
-- ^Schur factorization. (u,s) = schur m ==> m == u * s * u' where u is unitary and s is a Schur matrix. A complex Schur matrix is upper triangular. A real Schur matrix is upper triangular in 2x2 blocks.
schur m = let (u,s) = H.schur $ repa2hm m in (hm2repa u, hm2repa s)

schurS :: (Field t, Numeric t) => Array D DIM2 t -> (Array F DIM2 t, Array F DIM2 t)
schurS m = let (u,s) = H.schur $ repa2hmS m in (hm2repa u, hm2repa s)

schurSIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, Array F DIM2 t)
schurSIO m = do 
  (u,s) <- H.schur <$> repa2hmSIO m 
  return (hm2repa u, hm2repa s)

schurP :: (Field t, Numeric t, Monad m) => Array D DIM2 t -> m (Array F DIM2 t, Array F DIM2 t)
schurP m = do 
  (u,s) <- H.schur <$> repa2hmP m 
  return (hm2repa u, hm2repa s)

schurPIO :: (Field t, Numeric t) => Array D DIM2 t -> IO (Array F DIM2 t, Array F DIM2 t)
schurPIO m = do 
  (u,s) <- H.schur <$> repa2hmPIO m 
  return (hm2repa u, hm2repa s)
