{-# LANGUAGE
    DataKinds
  , ExistentialQuantification
  , FlexibleInstances
  , PolyKinds
  , RankNTypes
  , ScopedTypeVariables
  , TypeFamilies
  , TypeOperators
  , UndecidableInstances
  #-}

module Numeric.LinearAlgebra.Repa.Internal.Refined where

import Data.Array.Repa hiding (reshape)
import Data.Proxy
import Data.Reflection hiding (Z)
import Data.Coerce
import Data.Type.Equality
import GHC.TypeLits

newtype RArray refs t sh e = RArray (Array t sh e)
newtype IArray (ns :: [Nat]) refs t sh e = IArray (Array t sh e)

class ReifyShape sh where
  reifyShape :: sh -> (forall (ns :: [Nat]). Proxy ns -> a) -> a

instance ReifyShape Z where
  reifyShape Z f = f p
    where p :: Proxy '[] = Proxy

instance ReifyShape sh => ReifyShape (sh :. Int) where
  reifyShape (sh :. i) f = reifyNat (fromIntegral i) $ \p -> reifyShape sh $ \ps -> f (p `appendType` ps)

appendType :: proxy1 n -> proxy2 ns -> Proxy (n ': ns)
appendType _ _ = Proxy

setSize :: Proxy ns -> RArray r t sh e -> IArray ns r t sh e
setSize _ = coerce

forgetSize :: IArray ns r t sh e -> RArray r t sh e
forgetSize = coerce

(.+.) :: Proxy a -> Proxy b -> Proxy (a+b)
Proxy .+. Proxy = Proxy

(.-.) :: Proxy a -> Proxy b -> Proxy (a+b)
Proxy .-. Proxy = Proxy

(.*.) :: Proxy a -> Proxy b -> Proxy (a*b)
Proxy .*. Proxy = Proxy

(.=.) :: Proxy a -> Proxy b -> Proxy (a==b)
Proxy .=. Proxy = Proxy
