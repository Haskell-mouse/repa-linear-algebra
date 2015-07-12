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
  ) where

import Numeric.LinearAlgebra.Repa.Conversion

import Data.Array.Repa hiding (rank)
import Data.Array.Repa.Repr.ForeignPtr
import qualified Numeric.LinearAlgebra.HMatrix as H
import Numeric.LinearAlgebra.HMatrix (Numeric, Field, LSDiv, Normed, Product, Vector)

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
