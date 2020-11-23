import bpy
import os, sys, argparse

def main():	
	parser = argparse.ArgumentParser()
    parser.add_argument("--render", action='store_true')
    if '--' in sys.argv:
        args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    else:
        args = parser.parse_args('')

    # render
    bpy.context.scene.render.resolution_percentage = 50
    bpy.context.scene.cycles.samples = 60
    bpy.context.scene.render.engine = 'CYCLES'
    if args.render:
	    bpy.ops.render.render()
	    bpy.data.images['Render Result'].save_render(filepath = os.path.dirname(__file__)+'/out/01_out.png')

if __name__ == "__main__":
    main()