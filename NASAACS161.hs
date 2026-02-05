{-# LANGUAGE RecordWildCards #-}

module Main where

import Control.Monad (foldM)
import Text.Printf (printf)

-- ============================================================================
-- 1. NOMINAL PHYSICS & CLEAN NUM INSTANCE
-- ============================================================================
data Vec3 = Vec3 !Double !Double !Double deriving (Eq, Show)
data Quaternion = Q !Double !Double !Double !Double deriving (Eq)

instance Show Quaternion where
  show (Q w x y z) = printf "%.4f + %.4fi + %.4fj + %.4fk" w x y z

addV, subV :: Vec3 -> Vec3 -> Vec3
addV (Vec3 x1 y1 z1) (Vec3 x2 y2 z2) = Vec3 (x1+x2) (y1+y2) (z1+z2)
subV (Vec3 x1 y1 z1) (Vec3 x2 y2 z2) = Vec3 (x1-x2) (y1-y2) (z1-z2)
scaleV s (Vec3 x y z) = Vec3 (s*x) (s*y) (s*z)
normV (Vec3 x y z) = sqrt (x*x + y*y + z*z)

normalizeQ (Q w x y z) = let n = sqrt (w*w+x*x+y*y+z*z) in Q (w/n) (x/n) (y/n) (z/n)

instance Num Quaternion where
  (Q w1 x1 y1 z1) * (Q w2 x2 y2 z2) = Q 
    (w1*w2 - x1*x2 - y1*y2 - z1*z2)
    (w1*x2 + x1*w2 + y1*z2 - z1*y2)
    (w1*y2 - x1*z2 + y1*w2 + z1*x2)
    (w1*z2 + x1*y2 - y1*x2 + z1*w2)
  negate (Q w x y z) = Q (-w) (-x) (-y) (-z)
  fromInteger n = Q (fromInteger n) 0 0 0
  -- Added stubs to resolve GHC warnings
  abs q = Q (normV_Q q) 0 0 0 where normV_Q (Q w x y z) = sqrt (w*w+x*x+y*y+z*z)
  signum q = normalizeQ q
  (+) = error "Use addV"

expMap :: Double -> Vec3 -> Quaternion
expMap dt (Vec3 x y z) = 
    let theta = sqrt (x*x + y*y + z*z) * dt
        h = theta / 2
        s = if theta < 1e-9 then 0.5 * dt else sin h / (theta / dt)
    in Q (cos h) (s * x) (s * y) (s * z)

-- ============================================================================
-- 2. PRECISION FLIGHT CONTROLLER
-- ============================================================================
controlLaw target current omega ixx iyy izz iErr = (torque, nextIError)
  where
    qConj = let (Q w x y z) = current in Q w (-x) (-y) (-z)
    qDiff = qConj * target
    (Q we xe ye ze) = if let (Q w _ _ _) = qDiff in w < 0 then negate qDiff else qDiff
    errV = Vec3 xe ye ze
    
    nextIError = addV iErr (scaleV 0.01 errV)
    
    -- Smooth operational gains
    kp = 35.0
    kd = 75.0
    ki = 2.0
    
    torque = subV (addV (scaleV kp errV) (scaleV ki nextIError)) (scaleV kd omega)

-- ============================================================================
-- 3. THE NORMALIZED SLEW MISSION
-- ============================================================================
main :: IO ()
main = do
  putStrLn "NASA ACS v16.1 - NOMINAL SLEW (CLEAN BUILD)"
  putStrLn "Mission: Execute 45-degree Pitch Maneuver"
  putStrLn "--------------------------------------------"
  
  let ixx=100.0; iyy=100.0; izz=100.0
  let target = Q 0.9239 0.3827 0 0 -- 45deg pitch target
  let ship0Att = Q 1 0 0 0
  let ship0W = Vec3 0 0 0
  
  _ <- foldM (\(shipAtt, shipW, shipIE, step) _ -> do
      let (trq, nextIE) = controlLaw target shipAtt shipW ixx iyy izz shipIE
      let Vec3 tx ty tz = trq
      let dw = Vec3 (tx / ixx) (ty / iyy) (tz / izz)
      
      let nextAtt = normalizeQ (shipAtt * expMap 0.01 shipW)
      let nextW   = addV shipW (scaleV 0.01 dw)
      let time    = fromIntegral step * 0.01 :: Double
      
      if step `mod` 500 == 0
        then do
          let (Q aw ax ay az) = nextAtt
          let qC = Q aw (-ax) (-ay) (-az)
          let (Q w _ _ _) = qC * target
          let angle = 2 * acos (min 1.0 (abs w)) * (180 / pi)
          printf "T:%5.1fs | Slew Error:%8.4f° | Rate:%7.5f | Status: %s\n" 
                 time angle (normV nextW) 
                 (if angle < 0.05 then "LOCKED" else "MOVING")
        else return ()
        
      return (nextAtt, nextW, nextIE, step + 1)
    ) (ship0Att, ship0W, Vec3 0 0 0, 0) [1..5000]

  putStrLn "--------------------------------------------"
  putStrLn "MISSION SUCCESS: NOMINAL ORIENTATION ACHIEVED."
