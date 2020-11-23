import os, sys, math, argparse # python modules
import bpy, mathutils # blender modules

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import myutil # my functions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--render", action='store_true')
    if '--' in sys.argv:
        args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    else:
        args = parser.parse_args('')


    myutil.set_tex_environment(
    	bpy.context.scene.world,
    	os.path.dirname(os.path.abspath(__file__))+'/../test_inputs/epping_forest_01_1k.hdr')

    # rotate cube
    bpy.data.objects["Cube"].rotation_euler = _euler = mathutils.Euler((0.0, 0.0, math.radians(30.0)), 'XYZ')

    myutil.set_camera_ydirection()

    # render
    bpy.context.scene.render.resolution_percentage = 50
    bpy.context.scene.cycles.samples = 60
    bpy.context.scene.render.engine = 'CYCLES'
    if args.render:
        bpy.ops.render.render()
        bpy.data.images['Render Result'].save_render(filepath = os.path.dirname(__file__)+'/out/02_out.png')

if __name__ == "__main__":    
	main()   
