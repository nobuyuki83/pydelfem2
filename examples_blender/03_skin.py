import bpy
import os, sys, argparse

def main():		
    parser = argparse.ArgumentParser()
    parser.add_argument("--render", action='store_true')
    if '--' in sys.argv:
        args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    else:
        args = parser.parse_args('')

    # remove initial cube
    bpy.data.objects.remove(bpy.data.objects["Cube"])  

    # add floor
    bpy.ops.mesh.primitive_plane_add(size=20)

    # set initial light location
    bpy.data.objects["Light"].location = (-4,4,4) 

    material_skin = bpy.data.materials.new('Skin')
    material_skin.use_nodes = True	    
    material_skin.node_tree.nodes["Principled BSDF"].inputs[0].default_value = (0.8, 0.2, 0.0, 1)
    material_skin.node_tree.nodes["Principled BSDF"].inputs[1].default_value = 0.5
    material_skin.node_tree.nodes["Principled BSDF"].inputs[2].default_value = (1.0, 0.1, 0.1)
    material_skin.node_tree.nodes["Principled BSDF"].inputs[3].default_value = (1.0, 1.0, 1.0, 1)
    material_skin.node_tree.nodes["Principled BSDF"].inputs[5].default_value = 0.2

    bpy.ops.mesh.primitive_monkey_add(location=(0,0,1), rotation=(0,0,45))
    bpy.ops.object.shade_smooth()   
    bpy.context.object.data.materials.append(material_skin)

    bpy.data.cameras["Camera"].lens = 90
    bpy.data.objects["Camera"].location[2] = 6.0

    # render
    bpy.context.scene.render.resolution_percentage = 50
    bpy.context.scene.cycles.samples = 60
    bpy.context.scene.render.engine = 'CYCLES'    
    if args.render:
        bpy.ops.render.render()
        bpy.data.images['Render Result'].save_render(filepath = os.path.dirname(__file__)+'/out/03_out.png')

if __name__ == "__main__":
    main()