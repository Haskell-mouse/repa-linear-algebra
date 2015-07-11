{-# LANGUAGE 
    FlexibleContexts
  , FlexibleInstances
  , TypeFamilies
  , TypeSynonymInstances 
  #-}

module Numeric.LinearAlgebra.Repa.Conversion
  ( Container
  , Element
  , HShape(..)
  -- * Vector conversion utilities
  , hv2repa
  , repa2hv
  , repa2hvS
  , repa2hvSIO
  , repa2hvP
  , repa2hvPIO
  -- * Matrix conversion utilities
  , hm2repa
  , repa2hm
  , repa2hmS
  , repa2hmSIO
  , repa2hmP
  , repa2hmPIO
  ) where

import Data.Array.Repa
import Data.Array.Repa.Repr.ForeignPtr
import Foreign.Storable
import Foreign.ForeignPtr
import qualified Numeric.LinearAlgebra.HMatrix as H
import Numeric.LinearAlgebra.HMatrix (Container, Element, Numeric)
import qualified Data.Vector.Storable as V

-- |Shape-polymorphic conversion.
class HShape sh where
  type HType sh :: * -> *
  toRepa      :: Numeric t => HType sh t   -> Array F sh t
  fromRepa    :: Numeric t => Array F sh t -> HType sh t
  fromRepaS   :: Numeric t => Array D sh t -> HType sh t
  fromRepaSIO :: Numeric t => Array D sh t -> IO (HType sh t)
  fromRepaP   :: (Numeric t, Monad m) => Array D sh t -> m (HType sh t)
  fromRepaPIO :: Numeric t => Array D sh t -> IO (HType sh t)

instance HShape DIM1 where
  type HType DIM1 = H.Vector
  toRepa      = hv2repa
  fromRepa    = repa2hv
  fromRepaS   = repa2hvS
  fromRepaSIO = repa2hvSIO
  fromRepaP   = repa2hvP
  fromRepaPIO = repa2hvPIO

instance HShape DIM2 where
  type HType DIM2 = H.Matrix
  toRepa      = hm2repa
  fromRepa    = repa2hm
  fromRepaS   = repa2hmS
  fromRepaSIO = repa2hmSIO
  fromRepaP   = repa2hmP
  fromRepaPIO = repa2hmPIO

-- Vector conversion utilities

hv2repa :: Storable t => H.Vector t -> Array F DIM1 t
-- ^O(1). Convert a HMatrix Vector to a Repa Array.
hv2repa hv = fromForeignPtr (ix1 ln) ptr
  where (ptr, ln) = V.unsafeToForeignPtr0 hv

repa2hv :: Storable t => Array F DIM1 t -> H.Vector t
-- ^O(1). Convert a Repa Array to a HMatrix Vector.
repa2hv r = V.unsafeFromForeignPtr0 (toForeignPtr r) ln
  where ln = size $ extent r

repa2hvS :: Storable t => Array D DIM1 t -> H.Vector t
-- ^Convert a delayed Repa Array to a HMatrix Vector, evaluating it sequentially.
repa2hvS = repa2hv.computeS

repa2hvSIO :: Storable t => Array D DIM1 t -> IO (H.Vector t)
-- ^O(1). Convert a Repa Array to a HMatrix Vector sequentially inside the IO monad.
repa2hvSIO r = do
    ptr <- mallocForeignPtrArray ln
    computeIntoS ptr r 
    return $ V.unsafeFromForeignPtr0 ptr ln
  where ln = size $ extent r

repa2hvP :: (Storable t, Monad m) => Array D DIM1 t -> m (H.Vector t)
-- ^Convert a delayed Repa Array to a HMatrix Vector, evaluating it in parallel.
repa2hvP = fmap repa2hv.computeP

repa2hvPIO :: Storable t => Array D DIM1 t -> IO (H.Vector t)
-- ^O(1). Convert a Repa Array to a HMatrix Vector in parallel inside the IO monad.
repa2hvPIO r = do
    ptr <- mallocForeignPtrArray ln
    computeIntoP ptr r 
    return $ V.unsafeFromForeignPtr0 ptr ln
  where ln = size $ extent r

-- Matrix conversion utilities

hm2repa 
  :: ( Storable t
     , Container V.Vector t
     , Element t
     ) 
  => H.Matrix t -> Array F DIM2 t
-- ^O(1). Convert a HMatrix Matrix to a Repa Array.
hm2repa hm = fromForeignPtr (ix2 r c) ptr
  where (ptr, _) = V.unsafeToForeignPtr0 $ H.flatten hm 
        (r  , c) = H.size hm

repa2hm :: Storable t => Array F DIM2 t -> H.Matrix t
-- ^O(1). Convert a Repa Array to a HMatrix Matrix.
repa2hm r = H.reshape c $ V.unsafeFromForeignPtr0 (toForeignPtr r) ln
  where ln = size e
        (_:c:[]) = listOfShape e
        e = extent r

repa2hmS :: Storable t => Array D DIM2 t -> H.Matrix t
-- ^Convert a delayed Repa Array to a HMatrix Matrix, evaluating it sequentially.
repa2hmS = repa2hm . computeS

repa2hmSIO :: Storable t => Array D DIM2 t -> IO (H.Matrix t)
-- ^O(1). Convert a Repa Array to a HMatrix Matrix sequentially inside the IO monad.
repa2hmSIO r = do
    ptr <- mallocForeignPtrArray ln
    computeIntoS ptr r 
    return . H.reshape c  $ V.unsafeFromForeignPtr0 ptr ln
  where ln = size e
        (_:c:[]) = listOfShape e
        e = extent r

repa2hmP :: (Storable t, Monad m) => Array D DIM2 t -> m (H.Matrix t)
-- ^Convert a delayed Repa Array to a HMatrix Matrix, evaluating it in parallel.
repa2hmP = fmap repa2hm . computeP

repa2hmPIO :: Storable t => Array D DIM2 t -> IO (H.Matrix t)
-- ^O(1). Convert a Repa Array to a HMatrix Matrix in parallel inside the IO monad.
repa2hmPIO r = do
    ptr <- mallocForeignPtrArray ln
    computeIntoP ptr r 
    return . H.reshape c  $ V.unsafeFromForeignPtr0 ptr ln
  where ln = size e
        (_:c:[]) = listOfShape e
        e = extent r


