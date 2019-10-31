import ffmpeg
import numpy as np
import contextlib

# Assumes frames is a numpy ndarray, where frames.shape = (n, w, h) where n is
# the number of frames, and w/h is the width/height of each frame.
def save_video(frames, out_path, scale=1):
	with video_stream(frames[0].shape, out_path) as vid_out:
		for f in frames:
			vid_out.write(f)

@contextlib.contextmanager
def video_stream(dims, out_path, scale=1):
	vid_out = (
	    ffmpeg
	    .input('pipe:', format='rawvideo',
	    	pix_fmt='gray8', s='{}x{}'.format(*(dims)))
	    .filter('scale', s='{}x{}'.format(dims[0] * scale, dims[1] * scale), flags='neighbor')
	    .output(out_path, pix_fmt='yuv420p')
	    .overwrite_output()
	    .run_async(pipe_stdin=True, quiet=True)
	)
	vid_out.write = lambda frame : vid_out.stdin.write(frame.astype(np.uint8).tobytes())
	yield vid_out
	vid_out.stdin.close()
	vid_out.wait()

def example():
	frames = [np.random.choice([0, 100, 255], size=(100, 100)) for _ in range(100)]
	save_video(frames, 'test2.mp4')