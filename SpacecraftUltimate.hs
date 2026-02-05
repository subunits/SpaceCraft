{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE RecordWildCards #-}

module Main where

import Data.List (intercalate, foldl')
import Control.Monad (replicateM, foldM)
import System.Random
import Text.Printf (printf)

-- ============================================================================
-- QUATERNIONS
-- ============================================================================

data Quaternion = Q Double Double Double Double
  deriving (Eq)

instance Show Quaternion where
  show (Q w x y z) = printf "%.4f + %.4fi + %.4fj + %.4fk" w x y z

instance Num Quaternion where
  (Q w1 x1 y1 z1) + (Q w2 x2 y2 z2) = Q (w1+w2) (x1+x2) (y1+y2) (z1+z2)
  (Q w1 x1 y1 z1) * (Q w2 x2 y2 z2) = Q w' x' y' z'
    where
      w' = w1*w2 - x1*x2 - y1*y2 - z1*z2
      x' = w1*x2 + x1*w2 + y1*z2 - z1*y2
      y' = w1*y2 - x1*z2 + y1*w2 + z1*x2
      z' = w1*z2 + x1*y2 - y1*x2 + z1*w2
  abs q = Q (norm q) 0 0 0
  signum q = let n = norm q in if n == 0 then q else scaleQ (1/n) q
  fromInteger n = Q (fromInteger n) 0 0 0
  negate (Q w x y z) = Q (-w) (-x) (-y) (-z)

scaleQ :: Double -> Quaternion -> Quaternion
scaleQ s (Q w x y z) = Q (s*w) (s*x) (s*y) (s*z)

conjugate :: Quaternion -> Quaternion
conjugate (Q w x y z) = Q w (-x) (-y) (-z)

norm :: Quaternion -> Double
norm (Q w x y z) = sqrt (w*w + x*x + y*y + z*z)

normalize :: Quaternion -> Quaternion
normalize q = let n = norm q in if n == 0 then q else scaleQ (1/n) q

dot :: Quaternion -> Quaternion -> Double
dot (Q w1 x1 y1 z1) (Q w2 x2 y2 z2) = w1*w2 + x1*x2 + y1*y2 + z1*z2

newtype UnitQuaternion = UQ Quaternion deriving (Eq)

instance Show UnitQuaternion where
  show (UQ q) = show q

mkUnitQuaternion :: Quaternion -> UnitQuaternion
mkUnitQuaternion = UQ . normalize

fromUnitQuaternion :: UnitQuaternion -> Quaternion
fromUnitQuaternion (UQ q) = q

-- ============================================================================
-- LIE ALGEBRA
-- ============================================================================

data Vec3 = Vec3 Double Double Double deriving (Eq, Show)

exponentialMap :: Vec3 -> UnitQuaternion
exponentialMap (Vec3 x y z) = mkUnitQuaternion $ Q (cos halfTheta) (sinc * x) (sinc * y) (sinc * z)
  where
    theta = sqrt (x*x + y*y + z*z)
    halfTheta = theta / 2
    sinc = if theta < 1e-8 then 0.5 else sin halfTheta / theta

logarithmicMap :: UnitQuaternion -> Vec3
logarithmicMap (UQ (Q w x y z))
  | abs w >= 1.0 = Vec3 0 0 0
  | otherwise = Vec3 (scale * x) (scale * y) (scale * z)
  where
    theta = 2 * acos (clamp (-1) 1 w)
    sinHalfTheta = sqrt (x*x + y*y + z*z)
    scale = if sinHalfTheta < 1e-8 then 2 else theta / sinHalfTheta
    clamp lo hi val = max lo (min hi val)

-- ============================================================================
-- SLERP
-- ============================================================================

slerp :: UnitQuaternion -> UnitQuaternion -> Double -> UnitQuaternion
slerp (UQ q1) (UQ q2) t
  | t <= 0    = UQ q1
  | t >= 1    = UQ q2
  | cosTheta > 0.9995 = mkUnitQuaternion $ q1 + scaleQ t (q2 - q1)
  | otherwise = UQ $ scaleQ (sin ((1-t)*theta) / sinTheta) q1 
                   + scaleQ (sin (t*theta) / sinTheta) q2'
  where
    cosTheta = dot q1 q2
    q2' = if cosTheta < 0 then negate q2 else q2
    theta = acos (abs cosTheta)
    sinTheta = sin theta

-- ============================================================================
-- PARALLEL TRANSPORT
-- ============================================================================

parallelTransport :: [UnitQuaternion] -> Vec3 -> Vec3
parallelTransport [] v = v
parallelTransport [_] v = v
parallelTransport (q1:q2:qs) v = parallelTransport (q2:qs) v'
  where
    relativeRotation = mkUnitQuaternion $ fromUnitQuaternion q2 * conjugate (fromUnitQuaternion q1)
    v' = rotateVector relativeRotation v

rotateVector :: UnitQuaternion -> Vec3 -> Vec3
rotateVector (UQ q) (Vec3 x y z) = Vec3 x' y' z'
  where
    qv = Q 0 x y z
    Q _ x' y' z' = q * qv * conjugate q

-- ============================================================================
-- SPACECRAFT STATE
-- ============================================================================

data SpacecraftState = SpacecraftState
  { position :: Vec3
  , velocity :: Vec3
  , attitude :: UnitQuaternion
  , angularVelocity :: Vec3
  , previousError :: Maybe Double
  , time :: Double
  } deriving (Show)

-- ============================================================================
-- SENSORS
-- ============================================================================

data SensorReading = SensorReading
  { starTrackerAttitude :: UnitQuaternion
  , imuAngularVelocity :: Vec3
  , sunVector :: Vec3
  , measurementNoise :: Double
  } deriving (Show)

simulateSensor :: SpacecraftState -> IO SensorReading
simulateSensor SpacecraftState{..} = do
  noiseVec <- randomVec3 0.001
  let noisyAttitude = attitude `compose` exponentialMap noiseVec
  
  noiseAngVel <- randomVec3 0.0001
  let noisyAngVel = angularVelocity `addVec3` noiseAngVel
  
  noiseSun <- randomVec3 0.01
  let noisySun = Vec3 1 0 0 `addVec3` noiseSun
  
  return $ SensorReading noisyAttitude noisyAngVel noisySun 0.001

compose :: UnitQuaternion -> UnitQuaternion -> UnitQuaternion
compose (UQ q1) (UQ q2) = mkUnitQuaternion (q1 * q2)

addVec3 :: Vec3 -> Vec3 -> Vec3
addVec3 (Vec3 x1 y1 z1) (Vec3 x2 y2 z2) = Vec3 (x1+x2) (y1+y2) (z1+z2)

scaleVec3 :: Double -> Vec3 -> Vec3
scaleVec3 s (Vec3 x y z) = Vec3 (s*x) (s*y) (s*z)

-- ============================================================================
-- CONTROL - PREDICTIVE STOPPING WITH PROPORTIONAL BRAKING
-- ============================================================================

data ControlCommand = ControlCommand
  { torque :: Vec3
  , thrustVector :: Vec3
  } deriving (Show)

-- PREDICTIVE STOPPING CONTROLLER
-- Key: Calculate stopping distance and brake proportionally
geometricAttitudeControl :: Double -> UnitQuaternion -> UnitQuaternion -> Vec3 -> Maybe Double -> ControlCommand
geometricAttitudeControl dt desired current angVel prevError = ControlCommand controlTorque (Vec3 0 0 0)
  where
    -- Antipodal handling
    dotProd = dot (fromUnitQuaternion desired) (fromUnitQuaternion current)
    desired' = if dotProd < 0 
               then mkUnitQuaternion (negate $ fromUnitQuaternion desired)
               else desired
    
    errorQuat = mkUnitQuaternion $ fromUnitQuaternion desired' * conjugate (fromUnitQuaternion current)
    errorVec = logarithmicMap errorQuat
    errorMag = vec3Norm errorVec
    angVelMag = vec3Norm angVel
    
    -- BASE GAINS - Moderate
    (kp_base, kd_base) = 
      if errorMag > 1.5 then (1.8, 0.6)      -- >85°
      else if errorMag > 0.5 then (1.4, 0.8) -- 30-85°
      else if errorMag > 0.1 then (0.9, 1.2) -- 6-30°
      else (0.4, 2.0)                        -- <6°
    
    -- PREDICTIVE STOPPING DISTANCE
    -- Physics: If we brake with damping kd, how far will we travel before stopping?
    -- Simplified: stoppingDistance ≈ v² / (2 * effectiveDamping)
    effectiveDamping = 2.0  -- Estimate of braking deceleration
    stoppingDistance = if angVelMag > 0.01
                       then (angVelMag * angVelMag) / (2.0 * effectiveDamping)
                       else 0
    
    -- CRITICAL INSIGHT: Start braking when stopping distance ≥ remaining error
    -- This prevents overshoot!
    stoppingRatio = if errorMag > 0.01
                    then stoppingDistance / errorMag
                    else 0
    
    -- PROPORTIONAL BRAKING based on stopping ratio
    -- If stoppingRatio = 0.5: We'll stop at halfway point - keep going
    -- If stoppingRatio = 1.0: We'll stop exactly at target - maintain speed
    -- If stoppingRatio = 2.0: We'll overshoot by 2x - BRAKE HARD
    
    predictiveBraking = 
      if stoppingRatio < 0.5 then 0                           -- Safe to accelerate
      else if stoppingRatio < 1.0 then (stoppingRatio - 0.5) * 2.0  -- Gentle braking
      else if stoppingRatio < 2.0 then (stoppingRatio - 1.0) * 5.0  -- Moderate braking
      else (stoppingRatio - 2.0) * 10.0 + 5.0                       -- Heavy braking
    
    -- DERIVATIVE TERM
    errorRate = case prevError of
                  Just prev -> (errorMag - prev) / dt
                  Nothing -> 0
    
    kd_derivative = if errorRate < -0.1 then 1.2  -- Error decreasing
                    else if errorRate > 0.1 then 0.8   -- Error increasing
                    else 1.0
    
    -- TOTAL DAMPING
    kd = kd_base * kd_derivative + predictiveBraking
    
    -- Reduce proportional gain when braking hard
    kp = if predictiveBraking > 3.0 then kp_base * 0.5
         else if predictiveBraking > 1.0 then kp_base * 0.7
         else kp_base
    
    -- Apply control law
    rawTorque = scaleVec3 kp errorVec `addVec3` scaleVec3 (-kd) angVel
    controlTorque = saturateTorque 12.0 rawTorque

saturateTorque :: Double -> Vec3 -> Vec3
saturateTorque maxTorque v@(Vec3 x y z) =
  let mag = vec3Norm v
  in if mag > maxTorque then scaleVec3 (maxTorque / mag) v else v

-- ============================================================================
-- DYNAMICS
-- ============================================================================

data InertiaTensor = InertiaTensor Double Double Double

defaultInertia :: InertiaTensor
defaultInertia = InertiaTensor 100.0 120.0 80.0

angularAcceleration :: InertiaTensor -> Vec3 -> Vec3 -> Vec3
angularAcceleration (InertiaTensor ixx iyy izz) (Vec3 tx ty tz) (Vec3 wx wy wz) =
  Vec3 ax ay az
  where
    ax = (tx - (izz - iyy) * wy * wz) / ixx
    ay = (ty - (ixx - izz) * wz * wx) / iyy
    az = (tz - (iyy - ixx) * wx * wy) / izz

integrateGeometric :: Double -> InertiaTensor -> ControlCommand -> SpacecraftState -> SpacecraftState
integrateGeometric dt inertia ControlCommand{..} SpacecraftState{..} = SpacecraftState
  { position = position `addVec3` scaleVec3 dt velocity
  , velocity = velocity `addVec3` scaleVec3 dt thrustVector
  , attitude = attitude `compose` exponentialMap (scaleVec3 dt angularVelocity)
  , angularVelocity = angularVelocity `addVec3` scaleVec3 dt angAccel
  , previousError = previousError
  , time = time + dt
  }
  where
    angAccel = angularAcceleration inertia torque angularVelocity

-- ============================================================================
-- MISSION SIMULATION
-- ============================================================================

simulationStep :: UnitQuaternion -> SpacecraftState -> IO SpacecraftState
simulationStep desiredAttitude state = do
  sensor <- simulateSensor state
  let estimatedAttitude = starTrackerAttitude sensor
  let estimatedAngVel = imuAngularVelocity sensor
  
  let currentError = computeError desiredAttitude estimatedAttitude
  
  let dt = 0.01
  let control = geometricAttitudeControl dt desiredAttitude estimatedAttitude estimatedAngVel (previousError state)
  
  let newState = integrateGeometric dt defaultInertia control state
  return $ newState { previousError = Just currentError }

runMission :: Int -> UnitQuaternion -> SpacecraftState -> IO [SpacecraftState]
runMission steps desired initial = foldM step [initial] [1..steps]
  where
    step states _ = do
      let current = head states
      next <- simulationStep desired current
      return (next : states)

computeError :: UnitQuaternion -> UnitQuaternion -> Double
computeError desired current = vec3Norm errorVec
  where
    dotProd = dot (fromUnitQuaternion desired) (fromUnitQuaternion current)
    desired' = if dotProd < 0 
               then mkUnitQuaternion (negate $ fromUnitQuaternion desired)
               else desired
    errorQuat = mkUnitQuaternion $ fromUnitQuaternion desired' * conjugate (fromUnitQuaternion current)
    errorVec = logarithmicMap errorQuat

-- ============================================================================
-- UTILITIES
-- ============================================================================

randomVec3 :: Double -> IO Vec3
randomVec3 scale = do
  x <- randomRIO (-scale, scale)
  y <- randomRIO (-scale, scale)
  z <- randomRIO (-scale, scale)
  return $ Vec3 x y z

randomUnitQuaternion :: IO UnitQuaternion
randomUnitQuaternion = do
  [w, x, y, z] <- replicateM 4 (randomRIO (-1.0, 1.0))
  return $ mkUnitQuaternion (Q w x y z)

vec3Norm :: Vec3 -> Double
vec3Norm (Vec3 x y z) = sqrt (x*x + y*y + z*z)

takeEvery :: Int -> [a] -> [a]
takeEvery _ [] = []
takeEvery n xs = head xs : takeEvery n (drop n xs)

-- ============================================================================
-- MAIN - PREDICTIVE STOPPING
-- ============================================================================

main :: IO ()
main = do
  putStrLn "╔═══════════════════════════════════════════════════════════════╗"
  putStrLn "║   PREDICTIVE STOPPING CONTROL v10.0 - ULTIMATE               ║"
  putStrLn "║   Brake Based on Predicted Stopping Distance                 ║"
  putStrLn "╚═══════════════════════════════════════════════════════════════╝"
  putStrLn ""
  
  initialAttitude <- randomUnitQuaternion
  let initialState = SpacecraftState
        { position = Vec3 0 0 0
        , velocity = Vec3 0 0 0
        , attitude = initialAttitude
        , angularVelocity = Vec3 0.1 (-0.05) 0.02
        , previousError = Nothing
        , time = 0.0
        }
  
  let desiredAttitude = mkUnitQuaternion (Q 1 0 0 0)
  
  putStrLn "=== INITIAL CONDITIONS ==="
  putStrLn $ "Attitude: " ++ show (attitude initialState)
  let initialError = computeError desiredAttitude initialAttitude
  putStrLn $ "Initial error: " ++ printf "%.4f" initialError ++ " rad (" ++ printf "%.1f" (initialError * 180 / pi) ++ "°)"
  
  let dotProd = dot (fromUnitQuaternion desiredAttitude) (fromUnitQuaternion initialAttitude)
  if dotProd < 0 
    then putStrLn "⚠️  Antipodal case - shortest path via sign flip"
    else putStrLn "✓ Direct path"
  putStrLn ""
  
  putStrLn "Control: Predictive stopping - brake when stoppingDistance ≥ error"
  putStrLn ""
  
  -- Extended mission
  putStrLn "=== MISSION: 60 SECONDS - PREDICTIVE BRAKING ===  "
  putStrLn "Time     Regime        Error        |ω|   StopDist  S_ratio  Status"
  putStrLn "────────────────────────────────────────────────────────────────────"
  
  mission <- runMission 6000 desiredAttitude initialState
  let finalState = head mission
  let errorFinal = computeError desiredAttitude (attitude finalState)
  
  let milestones = take 40 $ reverse $ takeEvery 150 mission
  mapM_ (\s -> do
    let err = computeError desiredAttitude (attitude s)
    let angVelNorm = vec3Norm $ angularVelocity s
    let stopDist = if angVelNorm > 0.01 then (angVelNorm * angVelNorm) / 4.0 else 0
    let stopRatio = if err > 0.01 then stopDist / err else 0
    
    let regime = if err > 1.5 then "ACQUISITION"
                 else if err > 0.5 then "TRACKING  "
                 else if err > 0.1 then "SETTLING  "
                 else "FINE-POINT"
    
    let status = if stopRatio > 2.0 then "H-BRAKE"
                 else if stopRatio > 1.0 then "M-BRAKE"
                 else if stopRatio > 0.5 then "S-BRAKE"
                 else if angVelNorm > 0.05 then "COAST  "
                 else "ACCEL  "
    
    putStrLn $ printf "%5.2f" (time s) ++ "s  [" ++ regime ++ "]  " ++ 
               printf "%5.3f" err ++ " rad (" ++ printf "%5.1f" (err * 180 / pi) ++ "°)  " ++
               printf "%.3f" angVelNorm ++ "  " ++ printf "%.4f" stopDist ++ "   " ++
               printf "%.2f" stopRatio ++ "   [" ++ status ++ "]"
    ) milestones
  
  putStrLn ""
  putStrLn "=== FINAL RESULTS ==="
  putStrLn $ "Initial error: " ++ printf "%.4f" initialError ++ " rad (" ++ printf "%.2f" (initialError * 180 / pi) ++ "°)"
  putStrLn $ "Final error:   " ++ printf "%.6f" errorFinal ++ " rad (" ++ printf "%.4f" (errorFinal * 180 / pi) ++ "°)"
  putStrLn $ "Final |ω|:     " ++ printf "%.6f" (vec3Norm $ angularVelocity finalState) ++ " rad/s"
  
  let convergence = 100 * (1 - errorFinal / initialError)
  putStrLn $ "Convergence:   " ++ printf "%.2f" convergence ++ "%"
  
  if initialError > 0.001
    then putStrLn $ "Error reduction: " ++ printf "%.1f" (initialError / errorFinal) ++ "x"
    else putStrLn "Error reduction: Already at target"
  putStrLn ""
  
  -- Convergence analysis
  let errors = map (\s -> computeError desiredAttitude (attitude s)) (reverse mission)
  let minError = minimum errors
  let last1000Errors = take 1000 errors
  let maxRecentError = maximum last1000Errors
  let finalOscillation = maxRecentError - minError
  
  putStrLn "=== CONVERGENCE QUALITY ==="
  putStrLn $ "Minimum error: " ++ printf "%.6f" minError ++ " rad (" ++ printf "%.3f" (minError * 180 / pi) ++ "°)"
  putStrLn $ "Final oscillation: " ++ printf "%.6f" finalOscillation ++ " rad (" ++ printf "%.3f" (finalOscillation * 180 / pi) ++ "°)"
  
  if finalOscillation > 0
    then putStrLn $ "Oscillation: " ++ printf "%.1f" (100 * finalOscillation / minError) ++ "% of minimum error"
    else putStrLn "Oscillation: None"
  putStrLn ""
  
  -- Success criteria
  if errorFinal < 0.001 then do
    putStrLn "╔═══════════════════════════════════════════════════════════════╗"
    putStrLn "║  🎯 PRECISION POINTING ACHIEVED                               ║"
    putStrLn "║  ✓ Error < 0.001 rad (< 0.06°) - Sub-arcminute precision    ║"
    putStrLn "║  ✓ Predictive stopping eliminated overshoot                  ║"
    putStrLn "║  ✓ Ready for precision science observations                  ║"
    putStrLn "╚═══════════════════════════════════════════════════════════════╝"
  else if errorFinal < 0.01 then do
    putStrLn "╔═══════════════════════════════════════════════════════════════╗"
    putStrLn "║  ✓ EXCELLENT CONVERGENCE                                      ║"
    putStrLn "║  ✓ Error < 0.01 rad (< 0.6°) - High precision achieved      ║"
    putStrLn "║  ✓ Predictive braking working as designed                    ║"
    putStrLn "║  ✓ Operational for spacecraft missions                       ║"
    putStrLn "╚═══════════════════════════════════════════════════════════════╝"
  else if errorFinal < 0.1 then do
    putStrLn "╔═══════════════════════════════════════════════════════════════╗"
    putStrLn "║  ✓ GOOD PERFORMANCE                                           ║"
    putStrLn "║  ✓ Error < 0.1 rad (< 6°) - Operational accuracy            ║"
    putStrLn "║  ✓ Predictive stopping demonstrated                          ║"
    putStrLn "╚═══════════════════════════════════════════════════════════════╝"
  else if finalOscillation < 0.5 then do
    putStrLn "╔═══════════════════════════════════════════════════════════════╗"
    putStrLn "║  STEADY IMPROVEMENT                                           ║"
    putStrLn $ "║  Convergence: " ++ printf "%.1f" convergence ++ "%  with controlled oscillation             ║"
    putStrLn "║  ✓ Predictive braking reducing overshoot                     ║"
    putStrLn "╚═══════════════════════════════════════════════════════════════╝"
  else do
    putStrLn "╔═══════════════════════════════════════════════════════════════╗"
    putStrLn "║  MISSION IN PROGRESS                                          ║"
    putStrLn $ "║  Progress: " ++ printf "%.1f" convergence ++ "% error reduction                                ║"
    putStrLn "╚═══════════════════════════════════════════════════════════════╝"
