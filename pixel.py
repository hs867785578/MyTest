import gym
import matplotlib.pyplot as plt
from gym.wrappers.pixel_observation import PixelObservationWrapper
from gym.wrappers.gray_scale_observation import GrayScaleObservation

# raw_env = RawPendulumEnv(g=9.81, biased=True, damping=1.0)
env = gym.make('Pendulum-v1', g=9.81)
env.reset()
env = PixelObservationWrapper(env,pixels_only=True)
env = GrayScaleObservation(env['pixels'])
# print(state)
state = env.reset()
print(state)
plt.imshow(state)