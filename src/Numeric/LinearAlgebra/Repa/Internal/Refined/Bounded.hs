{-# LANGUAGE
    DataKinds
  , FlexibleInstances
  , PolyKinds
  , GADTs
  , GeneralizedNewtypeDeriving
  , ScopedTypeVariables
  , TypeFamilies
  , TypeOperators
  #-}

module Numeric.LinearAlgebra.Repa.Internal.Refined.Bounded where

import Prelude hiding (Bounded)

import Data.Array.Repa hiding (rank, (++))
import Data.Array.Repa.Shape 
import Data.Proxy
import GHC.Base (quotInt, remInt)
import GHC.TypeLits

stage :: String
stage = "Numeric.LinearAlgebra.Repa.Indexed.Bounded"

newtype Bounded (n :: Nat) = Bounded Int
  deriving (Eq)

instance KnownNat n => Show (Bounded n) where
  show (Bounded i) = show i ++ "/" ++ show (natVal p)
    where p :: Proxy n = Proxy

bounded :: KnownNat n => Proxy n -> Int -> Bounded n
bounded p n = if natVal p >= fromIntegral n then Bounded n else outOfBounds p 1 

outOfBounds :: KnownNat n => Proxy n -> Int -> a
outOfBounds p n = error $ stage ++ ".A bounded shape of dimension " ++ show (natVal p) ++ " cannot hold an array of size " ++ show n

instance (Shape sh, KnownNat n, 0 <= n) => Shape (sh :. Bounded n) where
  {-# INLINE [1] rank #-}
  rank (t :. _) = 1 + rank t
  {-# INLINE [1] zeroDim #-}
  zeroDim = zeroDim :. Bounded 0  
  {-# INLINE [1] unitDim #-}
  unitDim = unitDim :. bounded p 1
     where p :: Proxy n = Proxy
  {-# INLINE [1] intersectDim #-}
  intersectDim (t1 :. Bounded n1) (t2 :. Bounded n2) = intersectDim t1 t2 :. Bounded (min n1 n2)
  {-# INLINE [1] addDim #-}
  addDim (t1 :. Bounded n1) (t2 :. Bounded n2) = addDim t1 t2 :. bounded p (n1 + n2)
     where p :: Proxy n = Proxy
  {-# INLINE [1] size #-}
  size (t :. Bounded n) = n * size t
  {-# INLINE [1] sizeIsValid #-}
  sizeIsValid (t :. Bounded n)  
    | size t > 0 = n <= maxBound `div` size t
    | otherwise  = False
  {-# INLINE [1] toIndex #-}
  toIndex (ta :. Bounded a) (ti :. Bounded i) = toIndex ta ti * a + i
  {-# INLINE [1] fromIndex #-}
  fromIndex (t :. Bounded n) i = fromIndex t (i `quotInt` n) :. bounded p r
    where r | rank t == 0 = i
            | otherwise   = i `remInt` n
          p :: Proxy n = Proxy
  {-# INLINE [1] inShapeRange #-}
  inShapeRange (zs :. Bounded z) (t1 :. Bounded n1) (t2 :. Bounded n2)
     = (n2 >= z) && (n1 >= z) && inShapeRange zs t1 t2
  {-# NOINLINE listOfShape #-}
  listOfShape (t :. Bounded n) = n : listOfShape t
  {-# NOINLINE shapeOfList #-}
  shapeOfList ss = case ss of 
                     [] -> error $ stage ++ ".toList: empty list when converting to  (_ :. Int)"
                     (s:ss') -> shapeOfList ss' :. bounded p s
     where p :: Proxy n = Proxy

  {-# INLINE deepSeq #-}
  deepSeq (t :. Bounded n) x = deepSeq t (n `seq` x)
