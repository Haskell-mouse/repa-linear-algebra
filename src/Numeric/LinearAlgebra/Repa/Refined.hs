{-# LANGUAGE 
    DataKinds
  , NoMonomorphismRestriction
  , PolyKinds 
  , RankNTypes
  #-}

module Numeric.LinearAlgebra.Repa.Refined 
  ( RArray
  , ReifyShape
  , unknown
  , reifySize
  ) where

import Numeric.LinearAlgebra.Repa.Internal.Refined

import Data.Array.Repa
import Data.Coerce
import GHC.TypeLits

unknown :: Array t sh e -> RArray '[] t sh e
-- ^Turn an unrefined array into a refined array that feigns no hypotheses. 
unknown = coerce

reifySize 
  :: ( ReifyShape sh
     , Shape sh
     , Source t e
     ) 
  => RArray rs t sh e 
  -> (forall (ns :: [Nat]). IArray ns rs t sh e -> a) 
  -> a
reifySize rArr@(RArray arr) f = reifyShape (extent arr) $ \p -> f (setSize p rArr)





