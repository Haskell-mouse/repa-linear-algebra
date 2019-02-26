module Numeric.LinearAlgebra.Repa.Experimental where

import Data.Array.Repa hiding (rank)
import Data.Array.Repa.Repr.ForeignPtr
import Data.Bifunctor (first)
import Data.Proxy
import Numeric.LinearAlgebra.HMatrix
  (Complex, Field, LSDiv, Normed, Numeric, Product, RandDist(..), RealElement,
  Seed, Vector)
import qualified Numeric.LinearAlgebra.HMatrix as H
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)
import Prelude hiding (map)
import Unsafe.Coerce

import GHC.TypeLits

newtype Matrix r (m :: Nat) (n :: Nat) a  = Matrix{ _unMatrix :: Array r DIM2 a}

data SomeMatrix r a where
   SomeMatrix :: (KnownNat n, KnownNat m) => Matrix r m n a -> SomeMatrix r a

makeMatrix :: forall r a. Source r a => Array r DIM2 a -> SomeMatrix r a
makeMatrix array =
  let (Z :. y :. x) = extent array
      (Just k) = someNatVal $ fromIntegral y
      uk = case k of
        SomeNat (Proxy :: Proxy n) -> case (someNatVal $ fromIntegral y) of
          Just (SomeNat (Proxy :: Proxy m)) -> SomeMatrix $ (Matrix array :: Matrix r m n a)
  in undefined

withMatrix :: forall r a k . Source r a => Array r DIM2 a -> (forall x y . (KnownNat x, KnownNat y) => Matrix r x y a -> k) -> k
withMatrix array f =
  let (Z :. y :. x) = extent array
      (Just k) = someNatVal $ fromIntegral y
  in case k of
        SomeNat (Proxy :: Proxy n) -> case (someNatVal $ fromIntegral y) of
          Just (SomeNat (Proxy :: Proxy m)) -> f (Matrix array :: Matrix r m n a)

withMatrix2 :: forall r a k . Source r a => Array r DIM2 a -> Array r DIM2 a -> (forall x y x1 y1 . (KnownNat x, KnownNat y, KnownNat x1, KnownNat y1) => Matrix r x y a -> Matrix r x1 y1 a -> k) -> k
withMatrix2 array array2 f =
  let (Z :. y :. x) = extent array
      (Z :. y1 :. x1) = extent array2
      (Just k) = someNatVal $ fromIntegral y
      (Just m) = someNatVal $ fromIntegral x
      (Just l) = someNatVal $ fromIntegral x1
      (Just n) = someNatVal $ fromIntegral y1
  in case k of
        SomeNat (Proxy :: Proxy n) -> case m of
          SomeNat (Proxy :: Proxy m) -> case l of
            SomeNat (Proxy :: Proxy l) -> case n of
              SomeNat (Proxy :: Proxy z) -> f (Matrix array :: Matrix r m n a) (Matrix array :: Matrix r l z a)

type ResOfMul j h = Matrix

mulM :: forall x1 x2 y1 y2 . (KnownNat x1, KnownNat x2, KnownNat y1, KnownNat y2) => Matrix D x1 y1 Double -> Matrix D x2 y2 Double -> Matrix F x1 y2 Double
mulM mk1 mk2 =
  let res@(AForeignPtr (Z :. y :. x) _ _) = (_unMatrix mk1) `mulS` (_unMatrix mk2)
  in if (natVal (Proxy @x1) == fromIntegral x) && (natVal (Proxy @y1) == fromIntegral y)
     then Matrix $ res
     else (error "test")

mulF :: forall x1 x2 y1 y2 . (KnownNat x1, KnownNat x2, KnownNat y1, KnownNat y2) => Matrix F x1 y1 Double -> Matrix F x2 y2 Double -> Matrix F x1 y2 Double
mulF mk1 mk2 =
  let res@(AForeignPtr (Z :. y :. x) _ _) = (_unMatrix mk1) `mul` (_unMatrix mk2)
  in if (natVal (Proxy @x1) == fromIntegral x) && (natVal (Proxy @y1) == fromIntegral y)
     then Matrix $ res
     else (error "test")

transposeM :: forall x1 y1 . (KnownNat x1, KnownNat y1) => Matrix D x1 y1 Double -> Matrix D y1 x1 Double
transposeM m = Matrix $ transpose $ _unMatrix m

pFive :: forall x y . (KnownNat y, KnownNat x) => Matrix D x y Double -> Matrix D x y Double
pFive arr = Matrix $ map (+5.0) $ _unMatrix arr

sMult :: SomeMatrix D Double -> SomeMatrix D Double -> SomeMatrix F Double
sMult matr1 matr2 =
  case matr1 of
      (SomeMatrix mat1) -> case matr2 of
          (SomeMatrix mat2)  -> SomeMatrix $ mat1 `mulM` mat2

sMult2 :: Array D DIM2 Double -> Array D DIM2 Double -> Array F DIM2 Double
sMult2 arr1 arr2 = withMatrix2 arr1 arr2 (\a1 a2 -> _unMatrix $ mulM a1 a2)

mulk :: SomeMatrix D Double -> SomeMatrix D Double -> SomeMatrix F Double
mulk matr1 matr2 =
    case matr1 of
      (SomeMatrix mat1) -> case matr2 of
          (SomeMatrix mat2)  ->
              let mkm1 = transposeM mat1
                  mkm2 = mat2
              in SomeMatrix $ (mat1 `mulM` mat2) `mulF` (mkm1 `mulM` mkm2)
