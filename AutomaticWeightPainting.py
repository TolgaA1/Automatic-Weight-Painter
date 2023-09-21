import bpy
import bmesh
import mathutils
from mathutils import Vector
import math
from math import sqrt
from math import inf
from math import log
from math import acos
import time

#Welcome to the automatic weight painting Blender project
#- You will more than likely get a red question mark for the python file as blender does not update file path automatically. This means that
#the python file needs to be reloaded.

#- In the code section, to reload the file click on Text -> Open -> Find the directory of this project folder
#-> Open AutomaticWeightPainting.py under Python File folder.

#- Then click on the run button that looks like the play button in the code section. Will be at the top.

#- You will see a tab (already) open titled 'Tolga's Automatic Weight Paint System' in the left window that is the 3D viewport.

#- You can configure all the settings in the panel

#This system was built using Blender 3.1.0, although shouldn't be an issue, make sure to double check the version of Blender running.


def add_armature_modifier(mesh_object, selected_armature):
    #this function adds an armature modifier with selected armature to selected mesh
    #will replace already existing armature modifier if it exists
    is_armature_modifier_replaced = False 
    for modifier in mesh_object.modifiers:
        if modifier.type == "ARMATURE":
            mesh_object.modifiers.remove(mesh_object.modifiers.get(modifier.name))
            mesh_armature_modifier = mesh_object.modifiers.new('Armature', 'ARMATURE')
            mesh_armature_modifier.object = selected_armature
            is_armature_modifier_replaced = True
            
    if is_armature_modifier_replaced == False: 
            mesh_armature_modifier = mesh_object.modifiers.new('Armature', 'ARMATURE')
            mesh_armature_modifier.object = selected_armature

#function that resets all weighting    
def clean_weight_paint(vert):
    for group in vert.groups:
        group.weight = 0    
        
#function to get the distance between two points  
def get_distance_point(point, point2):
    distance = sqrt(
        (point[0] - point2[0]) * (point[0] - point2[0]) +
        (point[1] - point2[1]) * (point[1] - point2[1]) +
        (point[2] - point2[2]) * (point[2] - point2[2])
    )    
    return distance
      

#Automatic weight painting function for the spherical method
def auto_weight_paint_sphere_proximity_based(radius_offset, show_visual_shapes, radius_calc_method):
    time_start = time.time()
    #gets all the selected objects.
    selected_objects = [o for o in bpy.context.scene.objects if o.select_get()]
    print(radius_calc_method)     
    #Getting the selected armature, mesh and getting all the bones from the selected armature.
    #If the selected is not in intended order, switch mesh and armature around
    if selected_objects[0].type != 'MESH':
        selected_armature = selected_objects[0]
        mesh_object = bpy.data.objects[selected_objects[1].name]
        mesh_object.parent = bpy.data.objects[selected_objects[0].name]  
    else:
        selected_armature = selected_objects[1]
        mesh_object = bpy.data.objects[selected_objects[0].name]   
        mesh_object.parent = bpy.data.objects[selected_objects[1].name]  
    add_armature_modifier(mesh_object, selected_armature)
        
    bpy.context.view_layer.objects.active = selected_armature
  
    all_armature_bones = bpy.context.active_object.pose.bones
    bpy.ops.object.select_all(action='DESELECT')    

    #Delete all vertex groups to clear existing weights and groups
    mesh_object.vertex_groups.clear()
        
    for bone in all_armature_bones:
        #This will go through all the bones and will get their location. Will create a sphere around each bone
        #now will generate a radius for the influence sphere
        mesh_object_data = mesh_object.data
        weight_paint_radius =  0.3
        #List comphrehension to create a vertex list which will be looped through that are within 0.5 radius range of the bone.
        vert_list = [vert for vert in mesh_object_data.vertices if get_distance_point(mesh_object.matrix_world @ vert.co, bone.center) < 0.5]
        
        #Now it will get the nearest surface to current bone by getting minimum radius required
        #Check whether the radius calculation method was either nearest surface or average
        if radius_calc_method == "NearestSurface":
            for vert in vert_list:     
                #get distance between the bone center and vertex in this iteration.
                mesh_matrix = mesh_object.matrix_world
                vert_coord = vert.co
                #This gets the world coordinates of the vertex.                
                point = mesh_matrix @ vert_coord
                bone_coord = bone.center             
                distance_center = sqrt(
                    (point[0] - bone.center[0]) * (point[0] - bone.center[0]) +
                    (point[1] - bone.center[1]) * (point[1] - bone.center[1]) +
                    (point[2] - bone.center[2]) * (point[2] - bone.center[2])
                )
                if distance_center < weight_paint_radius:
                    weight_paint_radius = distance_center
                        
            weight_paint_radius += radius_offset + (-1/(bone.length+2))+0.5   

        elif radius_calc_method == "AverageRadius":
            distance_to_bone_values = []
            for vert in vert_list:
                #get distance between the bone center and vertex and append the distance of each vertex to the list
                mesh_matrix = mesh_object.matrix_world
                vert_coord = vert.co
                #This gets the world coordinates of the vertex.
                point = mesh_matrix @ vert_coord
                bone_coord = bone.center
                cur_distance = sqrt(
                    (point[0] - bone.center[0]) * (point[0] - bone.center[0]) +
                    (point[1] - bone.center[1]) * (point[1] - bone.center[1]) +
                    (point[2] - bone.center[2]) * (point[2] - bone.center[2])
                )
                if cur_distance < weight_paint_radius:
                    distance_to_bone_values.append(cur_distance)
            #Get average of all distance appended to the list and offset by the bone length
            if len(distance_to_bone_values) > 0:
                average_dist = sum(distance_to_bone_values)/len(distance_to_bone_values)
                weight_paint_radius = (average_dist)+ (-1/(bone.length+2))+0.5   
            else:
                weight_paint_radius = 0.1 + (bone.length*0.2)
                             
        #creating all visual influence spheres, one for tail, head and center bone if the user ticked this option in the UI
        if show_visual_shapes == True: 
            bpy.ops.mesh.primitive_uv_sphere_add(radius=weight_paint_radius, enter_editmode=False, align='WORLD', location= bone.tail, scale=(1, 1, 1))
            sphere_detection_mesh_tail = bpy.context.active_object
            
            bpy.ops.mesh.primitive_uv_sphere_add(radius=weight_paint_radius, enter_editmode=False, align='WORLD', location= bone.head, scale=(1, 1, 1))
            sphere_detection_mesh_head = bpy.context.active_object

                        
            bpy.ops.mesh.primitive_uv_sphere_add(radius=weight_paint_radius, enter_editmode=False, align='WORLD', location= bone.center, scale=(1, 1, 1))
            sphere_detection_mesh_center = bpy.context.active_object

        #Creates a new vertex group that belongs to the current bone.
        current_vertex_group = mesh_object.vertex_groups.new(name=bone.name)
            
        for vert in vert_list:
            #For each bone, go through all the vertices of the mesh and check whether the vertice is inside the influence sphere radius
            mesh_matrix = mesh_object.matrix_world
            vert_coord = vert.co
            point = mesh_matrix @ vert_coord
            distance_tail = get_distance_point(point, bone.tail)
            distance_head = get_distance_point(point, bone.head)
            distance_center = get_distance_point(point, bone.center)
            
            #if the vertex is in the influence sphere, assign it the vertex group of the bone with 100% influence
            if distance_center < weight_paint_radius or distance_tail < weight_paint_radius or distance_head < weight_paint_radius:
                current_vertex_group.add([vert.index], 1.0, 'REPLACE')
            
            
    time_end = time.time()
    print("Run time: ",time_end - time_start)

#-----------------------------------------------------------------------------------------------

def CylTest_CapsFirst(pt1, pt2, lengthsq, radius_sq, testpt ):
    #Code https://www.flipcode.com/archives/Fast_Point-In-Cylinder_Test.shtml from was used and adapted here
    dx = pt2.x - pt1.x
    dy = pt2.y - pt1.y    
    dz = pt2.z - pt1.z

    pdx = testpt.x - pt1.x		
    pdy = testpt.y - pt1.y
    pdz = testpt.z - pt1.z

    # Dot the d and pd vectors to see if point lies behind the 
    # cylinder cap at pt1.x, pt1.y, pt1.z

    dot = pdx * dx + pdy * dy + pdz * dz;

    # If dot is less than zero the point is behind the pt1 cap.
    # If greater than the cylinder axis line segment length squared
    # then the point is outside the other end cap at pt2.

    if dot < 0.0 or dot > lengthsq:
        return -1.0
    else:
        #If the point lies within the parallel caps, then find the distance squared from the point to line
        dsq = (pdx*pdx + pdy*pdy + pdz*pdz) - dot*dot/lengthsq
        #if the point is outside of the radius squared, return -1.0 which means it isnt inside the cylinder but if it is, then simply return
        #the squared value
        if( dsq > radius_sq ):
            return -1.0
        else:
            return dsq	
           

def cylinder_between(pos1, pos2, r):
  #Code from https://blender.stackexchange.com/questions/5898/how-can-i-create-a-cylinder-linking-two-points-with-python was used and adapted to this system
  #creates a cylinder from two given points and a radius
  dx = pos2.x - pos1.x
  dy = pos2.y - pos1.y
  dz = pos2.z - pos1.z   
  dist = math.sqrt(dx**2 + dy**2 + dz**2)

  #create the cylinder using Blender's function using the difference of xyz and halving it.    
  bpy.ops.mesh.primitive_cylinder_add(
      radius = r, 
      depth = dist,
      location = (dx/2 + pos1.x, dy/2 + pos1.y, dz/2 + pos1.z)   
  ) 
  
  #will get the rotation of of two points using atan2
  phi = math.atan2(dy, dx) 
  theta = math.acos(dz/dist) 
  #rotate the cylinder using phi and theta
  bpy.context.object.rotation_euler[1] = theta 
  bpy.context.object.rotation_euler[2] = phi 
  return bpy.context.object


#This is the function that handles the cylindrical method of the weight painting system.           
def auto_weight_paint_cylinder_proximity_based(radius_offset, show_visual_shapes, radius_calc_method,cylinder_length_offset):
    time_start = time.time() 
    #gets all the selected objects.
    selected_objects = [o for o in bpy.context.scene.objects if o.select_get()]
            
    #Getting the selected armature, mesh and getting all the bones from the selected armature. This is assuming the second picked object is the armature
    if selected_objects[0].type != 'MESH':
        selected_armature = selected_objects[0]
        mesh_object = bpy.data.objects[selected_objects[1].name]
        mesh_object.parent = bpy.data.objects[selected_objects[0].name]  
    else:
        selected_armature = selected_objects[1]
        mesh_object = bpy.data.objects[selected_objects[0].name]   
        mesh_object.parent = bpy.data.objects[selected_objects[1].name]  
    
    add_armature_modifier(mesh_object, selected_armature)
        
    bpy.context.view_layer.objects.active = selected_armature
    if cylinder_length_offset == 0:
        cylinder_length_offset = 1
    else:
        cylinder_length_offset = cylinder_length_offset * 2
        
    all_armature_bones = bpy.context.active_object.pose.bones
    bpy.ops.object.select_all(action='DESELECT')        
    #Delete all vertex groups to clear existing weights and groups
    mesh_object.vertex_groups.clear()
      
    for bone in all_armature_bones:
        #This will go through all the bones and will get their location. Will create a cylinder around each bone
        #now will generate a radius for the influence sphere
        mesh_object_data = mesh_object.data
        weight_paint_radius =  0.3
        bone_rotation = bone.matrix.to_euler()
        #Create a vertex list which will be looped through that are within a 0.5 radius of this bone
        vert_list = [vert for vert in mesh_object_data.vertices if get_distance_point(mesh_object.matrix_world @ vert.co, bone.center) < 0.5]
        #Now it will get the nearest surface to current bone by getting minimum radius required
        #Check whether the radius calculation method was either nearest surface or average
        if radius_calc_method == "NearestSurface":
            for vert in vert_list:     
                mesh_matrix = mesh_object.matrix_world
                vert_coord = vert.co
                #This gets the world coordinates of the vertex.
                point = mesh_matrix @ vert_coord
                bone_coord = bone.center             
                distance_center = sqrt(
                    (point[0] - bone.center[0]) * (point[0] - bone.center[0]) +
                    (point[1] - bone.center[1]) * (point[1] - bone.center[1]) +
                    (point[2] - bone.center[2]) * (point[2] - bone.center[2])
                )
                if distance_center < weight_paint_radius:
                    weight_paint_radius = distance_center
                        
            weight_paint_radius += radius_offset+ (-1/(bone.length+2))+0.5     
        
        elif radius_calc_method == "AverageRadius":
            distance_to_bone_values = []
            for vert in vert_list:
                #get distance between the bone center and vertex and append the distance of each vertex to the list
                mesh_matrix = mesh_object.matrix_world
                vert_coord = vert.co
                #This gets the world coordinates of the vertex.
                point = mesh_matrix @ vert_coord
                bone_coord = bone.center
                cur_distance = sqrt(
                    (point[0] - bone.center[0]) * (point[0] - bone.center[0]) +
                    (point[1] - bone.center[1]) * (point[1] - bone.center[1]) +
                    (point[2] - bone.center[2]) * (point[2] - bone.center[2])
                )
                if cur_distance < weight_paint_radius:
                    distance_to_bone_values.append(cur_distance)
            #now will get the average of all the values in this list
            if len(distance_to_bone_values) > 0:
                average_dist = sum(distance_to_bone_values)/len(distance_to_bone_values)
                weight_paint_radius = (average_dist)+ (-1/(bone.length+2))+0.5 
            else:
                weight_paint_radius = 0.1 + (bone.length*0.2)

        #calculating the vector direction of bone head and tail so the system can offset it so that the cylinder has larger length than the bone
        bone_direction_vector = (bone.tail - bone.head)/((0.05+ (-1/(bone.length+0.2))+5)/cylinder_length_offset)
        bone_tail_offset = bone.tail+bone_direction_vector
        new_bone_length = (bone_tail_offset - bone.head).length
        
        #create the visual cylinder shapes if it is selected in the user interface
        #using the function that creates a cylinder from two points.
        if show_visual_shapes == True:
            cylinder_detection_mesh = cylinder_between(bone_tail_offset,bone.head,weight_paint_radius)
            cylinder_detection_mesh.modifiers.new('Wireframe', 'WIREFRAME')

            
        current_vertex_group = mesh_object.vertex_groups.new(name=bone.name)
            
        for vert in vert_list:
            #For each bone, go through all the vertices of the mesh and check whether the vertice is inside the influence cylinder
            mesh_matrix = mesh_object.matrix_world
            vert_coord = vert.co
            point = mesh_matrix @ vert_coord
            bone_coord = bone.center
                
            distance_center = sqrt(
                (point[0] - bone.center.x) * (point[0] - bone.center.x) +
                (point[1] - bone.center.y) * (point[1] - bone.center.y) +
                (point[2] - bone.center.z) * (point[2] - bone.center.z)
            )  
            #if the point is inside the cylinder, give it a weighting of full value
            #check if the function result is not -1.0 since this means that the point is within the cylinder radius
            if CylTest_CapsFirst(bone_tail_offset, bone.head, new_bone_length*new_bone_length, weight_paint_radius*weight_paint_radius, point) != -1.0:
                current_vertex_group.add([vert.index], 1.0, 'REPLACE')
                
                
    time_end = time.time()
    print("Run time: ",time_end - time_start) 
    
#-----------------------------------------------------------------
            

#UI CODE FOR THE SCRIPT

bl_info = {
    "name": "Tolga's Automatic Weight Paint System",
    "author": "Tolga Arduc",
    "version": (0, 0, 1),
    "blender": (3, 1, 0),
    "location": "3D Viewport > Sidebar > Proximity Based Auto Weight Paint",
    "description": "An alternate automatic weight painting system that is more consistent than Blender's automatic weight painting system.",
    "category": "Development",
}


class AllProperties(bpy.types.PropertyGroup):
    #All properties for the UI get defined in this class.
    
    radius_offset_sphere: bpy.props.FloatProperty(
        name = "Sphere Radius Offset",
        default = 0.05,
        min = -0.25,
        max = 1,
        description = "Offset of proximity sphere radius"        
    )
    
    radius_offset_cylinder: bpy.props.FloatProperty(
        name = "Cylinder Radius Offset",
        default = 0.05,
        min = -0.25,
        max = 1,
        description = "The value of the offset for proximity cylinder radius"        
    )
    show_visual_shapes: bpy.props.BoolProperty(
        name="Show Influence Shapes",
        description="Will visually show areas influenced.",
        default = False
    )
    show_visual_cylinder: bpy.props.BoolProperty(
        name="Show Influence Shapes",
        description="Will visually show areas influenced.",
        default = False
    )
    show_visual_shapes_cylinder: bpy.props.BoolProperty(
        name="Show Influence Shapes",
        description="Will visually show areas influenced in the form of a cylinder.",
        default = False
    )
    
    radius_calc_method: bpy.props.EnumProperty(
        name="Radius Calculation Method",
        description="Choose how the radius is calculated.",
        items = [("NearestSurface", "Nearest Surface", ""),
                 ("AverageRadius", "Average Radius", "")
        ]
    )
    radius_calc_method_cylinder: bpy.props.EnumProperty(
        name="Radius Calculation Method",
        description="Choose how the radius is calculated.",
        items = [("NearestSurface", "Nearest Surface", ""),
                 ("AverageRadius", "Average Radius", "")
        ]
    ) 
    length_offset_cylinder: bpy.props.FloatProperty(
        name = "Cylinder Length Offset",
        default = 0.00,
        min = -5,
        max = 5,
        description = "The value of the offset for cylinder length"        
    )
    
#Defines the spherical method operator and the parameters it takes in.
class OBJECT_OT_auto_weight_proximity_sphere(bpy.types.Operator):
    """Testing"""
    bl_idname = "object.auto_weight_proximity_sphere"
    bl_label = "Automatic weight paint sphere"
    
    radius_offset_sphere: bpy.props.FloatProperty(
        name = "Sphere Radius Offset",
        default = 0.05,
        description = "Offset of proximity sphere's radius"
    )
    
    show_visual_shapes: bpy.props.BoolProperty(
        name="Show Influence Shapes",
        description="Will visually show areas influenced.",
        default = False
    )
    
    radius_calc_method: bpy.props.EnumProperty(
        name="Radius Calculation Method",
        description="Choose how the radius is calculated.",
        items = [("NearestSurface", "Nearest Surface", ""),
                 ("AverageRadius", "Average Radius", "")
        ]
    )
    #defining that the function for it will run whenever the operator is executed.
    def execute(self, context):
        auto_weight_paint_sphere_proximity_based(self.radius_offset_sphere, self.show_visual_shapes, self.radius_calc_method)
        return{"FINISHED"}

##Defines the cylindrical method operator and the parameters it takes in.    
class OBJECT_OT_auto_weight_proximity_cylinder(bpy.types.Operator):
    """Testing"""
    bl_idname = "object.auto_weight_proximity_cylinder"
    bl_label = "Automatic weight paint cylinder"
    
    radius_offset_cylinder: bpy.props.FloatProperty(
        name = "Cylinder Radius Offset",
        default = 0.05,
        description = "Offset of proximity cylinder's radius"
    )
    
    show_visual_shapes: bpy.props.BoolProperty(
        name="Show Influence Shapes",
        description="Will visually show areas influenced.",
        default = False
    )

    radius_calc_method: bpy.props.EnumProperty(
        name="Radius Calculation Method",
        description="Choose how the radius is calculated.",
        items = [("NearestSurface", "Nearest Surface", ""),
                 ("AverageRadius", "Average Radius", "")
        ]
    )
    
    length_offset_cylinder: bpy.props.FloatProperty(
        name = "Cylinder Length Offset",
        default = 0.00,
        min = -2,
        max = 2,
        description = "The value of the offset for cylinder length"        
    )
    
    def execute(self, context):
        auto_weight_paint_cylinder_proximity_based(self.radius_offset_cylinder, self.show_visual_shapes, self.radius_calc_method,self.length_offset_cylinder)
        return{"FINISHED"}

class VIEW3D_PT_my_custom_panel(bpy.types.Panel):

    # where to add the panel in the UI
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"  

    bl_category = "Auto Weight Paint"  # found in the Sidebar
    bl_label = "Tolga's Automatic Weight Paint System"  # found at the top of the Panel

    def draw(self, context):
        """define the layout of the panel"""
        row = self.layout.row()
        #defining the buttons and properties for the spherical method
        self.layout.prop(context.scene.properties, "show_visual_shapes")
        self.layout.prop(context.scene.properties, "radius_calc_method")
        row.label(text = "Sphere Based Proximity Method", icon ="MESH_UVSPHERE")
        self.layout.prop(context.scene.properties, "radius_offset_sphere")
        props = self.layout.operator("object.auto_weight_proximity_sphere", icon = "OUTLINER_OB_ARMATURE" , text="Automatic Sphere Weight Paint")
        props.radius_offset_sphere = context.scene.properties.radius_offset_sphere
        props.show_visual_shapes = context.scene.properties.show_visual_shapes
        props.radius_calc_method = context.scene.properties.radius_calc_method
        
        #Created a new row for the cylinderical method
        row = self.layout.row()
        #defining the buttons and properties for the spherical method
        row.label(text = "Cylinder Based Proximity Method", icon = "MESH_CYLINDER")
        self.layout.prop(context.scene.properties, "show_visual_shapes_cylinder")
        self.layout.prop(context.scene.properties, "radius_calc_method_cylinder")
        self.layout.prop(context.scene.properties, "radius_offset_cylinder")
        self.layout.prop(context.scene.properties, "length_offset_cylinder")
        props = self.layout.operator("object.auto_weight_proximity_cylinder", icon = "OUTLINER_OB_ARMATURE", text="Automatic Cylinder Weight Paint")
        props.radius_offset_cylinder = context.scene.properties.radius_offset_cylinder
        props.show_visual_shapes = context.scene.properties.show_visual_shapes_cylinder
        props.radius_calc_method = context.scene.properties.radius_calc_method_cylinder
        props.length_offset_cylinder = context.scene.properties.length_offset_cylinder

def register():
    bpy.utils.register_class(VIEW3D_PT_my_custom_panel)
    bpy.utils.register_class(OBJECT_OT_auto_weight_proximity_sphere)
    bpy.utils.register_class(OBJECT_OT_auto_weight_proximity_cylinder)
    bpy.utils.register_class(AllProperties)
    bpy.types.Scene.properties = bpy.props.PointerProperty(type = AllProperties)



def unregister():
    bpy.utils.unregister_class(VIEW3D_PT_my_custom_panel)
    bpy.utils.unregister_class(OBJECT_OT_auto_weight_proximity_sphere)
    bpy.utils.unregister_class(OBJECT_OT_auto_weight_proximity_cylinder)


if __name__ == "__main__":
    register()    
            
