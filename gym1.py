import gym
import torch
import numpy as np

env = gym.make("CartPole-v0")

state_dim = env.observation_space.shape[0]
# action_dim = env.action_space.shape[0]
# max_action = float(env.action_space.high[0])
print(env.action_space.n)
print(env.observation_space)


env2 = gym.make("Pendulum-v1")
print(env2.action_space)
print(env2.observation_space)
# print(action_dim)
# print(max_action)