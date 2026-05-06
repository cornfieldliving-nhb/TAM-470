import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
from pend_funcs import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
from pend_funcs import *

def pend_plot_gen_animation(t,y,animation_length=10, filename="pend.mp4"):
    
    filepath = os.path.join(os.getcwd(), filename)
    if os.path.exists(filepath):
        os.remove(filepath)
        
    assert y.shape[0] == 8, "Length of 'y' must be 8"
    assert t.shape[0] == y.shape[1], "Length of 't' must match the second dimension of 'y'"

    N = t.shape[0]
    fps = N / animation_length
    interval = 1000 / fps  # Calculate the correct interval in milliseconds

    fig, ax = plt.subplots()
    ax.set_xlim([-3, 3])
    ax.set_ylim([-3.5, 3.5])

    scatter_origin = ax.scatter(0, 0, c='blue', marker='.')
    scatter_p1 = ax.scatter([], [], c='blue', marker='o')
    scatter_p2 = ax.scatter([], [], c='blue', marker='o')
    line_p1 = ax.plot([], [], color='black')[0]
    line_p2 = ax.plot([], [], color='black')[0]
    trail = ax.scatter([], [], c='orange', marker='.', s=1)

    def animate(frame):
        
        scatter_p1.set_offsets(y[:2, frame])
        scatter_p2.set_offsets(y[2:4, frame])
        line_p1.set_data([0, y[0, frame]], [0, y[1, frame]])
        line_p2.set_data([y[2, frame], y[0, frame]], [y[3, frame], y[1, frame]])
        trail.set_offsets(y[2:4, :frame].T)
        

        if frame % 1000 == 0 and frame != 0:
            print(f"Processed {frame} frames")
        return scatter_p1, scatter_p2, line_p1, line_p2, trail

    anim = animation.FuncAnimation(fig, animate, frames=N, interval=interval, blit=True)
    anim.save(filename, writer='ffmpeg', fps=fps, dpi=150)
    plt.close(fig)

def pend_plot_gen_animation_fast(t, y, animation_length=10, fps=60, filename="pend_fast.mp4"):
    filepath = os.path.join(os.getcwd(), filename)
    if os.path.exists(filepath):
        os.remove(filepath)
    
    assert y.shape[0] == 8, "Length of 'y' must be 8"
    assert t.shape[0] == y.shape[1], "Length of 't' must match the second dimension of 'y'"

    N = t.shape[0]
    total_frames = int(animation_length * fps)
    step_size = max(1, N // total_frames)
    selected_frames = range(0, N, step_size)
    print(f"Selected {len(selected_frames)} frames out of {N}")

    fig, ax = plt.subplots()
    ax.set_xlim([-3, 3])
    ax.set_ylim([-3.5, 3.5])

    scatter_origin = ax.scatter(0, 0, c='blue', marker='.')
    scatter_p1 = ax.scatter([], [], c='blue', marker='o')
    scatter_p2 = ax.scatter([], [], c='blue', marker='o')
    line_p1, = ax.plot([], [], color='black')
    line_p2, = ax.plot([], [], color='black')
    trail = ax.scatter([], [], c='orange', marker='.', s=1)

    def animate(i):
        frame = selected_frames[i]
        scatter_p1.set_offsets(y[:2, frame])
        scatter_p2.set_offsets(y[2:4, frame])
        line_p1.set_data([0, y[0, frame]], [0, y[1, frame]])
        line_p2.set_data([y[2, frame], y[0, frame]], [y[3, frame], y[1, frame]])
        trail.set_offsets(y[2:4, :frame:step_size].T)

        if i % 100 == 0 and i != 0:
            print(f"Processed {i} frames")
        return scatter_p1, scatter_p2, line_p1, line_p2, trail

    anim = animation.FuncAnimation(fig, animate, frames=len(selected_frames), interval=1000/fps, blit=True)
    anim.save(filename, writer='ffmpeg', fps=fps, dpi=150)
    plt.close(fig)

def pend_plot_gen_animation_first_n_frames(t, y, n, fps=60, filename="pend_first_n.mp4"):
    filepath = os.path.join(os.getcwd(), filename)
    if os.path.exists(filepath):
        os.remove(filepath)
    
    assert t.shape[0] >= n, "Length of 't' must be greater than or equal to'n'"
    assert y.shape[0] == 8, "Length of 'y' must be 8"
    assert t.shape[0] == y.shape[1], "Length of 't' must match the second dimension of 'y'"

    N = t.shape[0]
    selected_frames = range(0, n)
    print(f"Selected first {len(selected_frames)} frames out of {N}")

    fig, ax = plt.subplots()
    ax.set_xlim([-3, 3])
    ax.set_ylim([-3.5, 3.5])

    scatter_origin = ax.scatter(0, 0, c='blue', marker='.')
    scatter_p1 = ax.scatter([], [], c='blue', marker='o')
    scatter_p2 = ax.scatter([], [], c='blue', marker='o')
    line_p1, = ax.plot([], [], color='black')
    line_p2, = ax.plot([], [], color='black')
    trail = ax.scatter([], [], c='orange', marker='.', s=1)

    def animate(i):
        frame = selected_frames[i]
        scatter_p1.set_offsets(y[:2, frame])
        scatter_p2.set_offsets(y[2:4, frame])
        line_p1.set_data([0, y[0, frame]], [0, y[1, frame]])
        line_p2.set_data([y[2, frame], y[0, frame]], [y[3, frame], y[1, frame]])
        trail.set_offsets(y[2:4, :frame].T)

        if i % 100 == 0 and i != 0:
            print(f"Processed {i} frames")
        return scatter_p1, scatter_p2, line_p1, line_p2, trail

    anim = animation.FuncAnimation(fig, animate, frames=len(selected_frames), interval=1000/fps, blit=True)
    anim.save(filename, writer='ffmpeg', fps=fps, dpi=150)
    plt.close(fig)

         
def pend_plot_final_instant(t,y,save_image=False):
    time = t.shape[0] - 1
    plt.clf()
    plt.xlim([-3,3])
    plt.ylim([-3.5,3.5])
    plt.scatter(0,0,c='blue',marker='.')
    plt.scatter(y[0,time],y[1,time],c='blue',marker='o')
    plt.plot([0,y[0,time]],[0,y[1,time]],color='black')
    plt.scatter(y[2,time],y[3,time],c='blue',marker='o')
    plt.plot([y[2,time],y[0,time]],[y[3,time],y[1,time]],color='black')
    plt.title("t = " + str(t[time]))
    if save_image:
        plt.savefig("final_instant.png")
    plt.show()