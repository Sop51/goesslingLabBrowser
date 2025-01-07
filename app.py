from flask import Flask, render_template, request, redirect, url_for, flash, session, jsonify
from forms import SubClusterForm, GeneForm, GroupForm, UploadFileForm, LoginForm, AnnotationForm, GroupSubplotForm, GeneSubplotForm, TimepointForm
import os
from dotenv import load_dotenv
import base64
import re
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import tempfile

# load environment variables from .env file
load_dotenv()

# initialize the app
app = Flask(__name__)
# configure the secret key
app.config['SECRET_KEY'] = os.getenv('SECRET_KEY', 'default_secret_key')
app.config['WTF_CSRF_ENABLED'] = False

# credentials from environment variables
USERNAME = os.getenv('USERNAME')
PASSWORD = os.getenv('PASSWORD')

# define and configure the upload folder destination
UPLOAD_FOLDER = 'uploads/'
app.config['UPLOAD_FILES_DEST'] = UPLOAD_FOLDER

# define and configure the subcluster folder destination
DATA_FOLDER = 'subclusterData/'
app.config['DATA_FILES_DEST'] = DATA_FOLDER

# create upload folder if it doesn't exist
if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)

# create the subcluster folder if it doesnt exist
if not os.path.exists(DATA_FOLDER):
    os.makedirs(DATA_FOLDER)

# create global variables to hold the loaded seurat objects
seurat_obj = None
Hepatocyte = None
Biliary_Epithelial_Cell = None
Endothelial_Cell = None
Hepatic_Stellate_Cell = None
Fibroblast = None
Macrophage = None
Monocyte = None
Erythrocyte = None
Neuron = None

# define a function to plot the gene UMAP
def plot_gene_umap(gene):
    pandas2ri.activate()
    # load the required libraries in R
    ro.r('library(Seurat)')
    ro.r('library(SeuratDisk)')
    # load in the seurat obj
    global seurat_obj
    # check if the gene exists in the seurat object
    gene_check = ro.r(f'"{gene}" %in% rownames(seurat_obj)')
    if not gene_check[0]:
        return None # return nothing if no
    # create the UMAP plot for the gene
    ro.r(f'gene_plot <- FeaturePlot(seurat_obj, features = "{gene}", order=TRUE, min.cutoff = "q10", repel = TRUE)')
    # create a temporary file to save the plot
    with tempfile.NamedTemporaryFile(delete=False, suffix='.png') as temp_file:
        temp_filename = temp_file.name
        # use the file name to create the plot
        ro.r(f'png("{temp_filename}", width = 800, height = 800, res = 150)')
        # print the plot to the device
        ro.r('print(gene_plot)')
        # finalize the plot
        ro.r('dev.off()')
    # read the plot file and convert to base64
    with open(temp_filename, "rb") as img_file:
        img_data = img_file.read()
        plot_url = base64.b64encode(img_data).decode('utf-8')
    # remove the temporary file after processing
    os.remove(temp_filename)
    # return the plot as a base64-encoded string to embed in HTML
    return f"data:image/png;base64,{plot_url}"

# define a function to plot the group umap
def plot_group_umap(group):
    pandas2ri.activate()
    # load the libraries in r
    ro.r('library(Seurat)')
    ro.r('library(SeuratDisk)')
    ro.r('library(ggplot2)')
    # load in the seurat obj
    global seurat_obj
    # check if the group exists in the seurat object metadata
    group_check = ro.r(f'"{group}" %in% colnames(seurat_obj@meta.data)')
    if not group_check[0]:
        return None
    
    # fix the error where the viewport is too small
    try:
        # create the UMAP plot for the group
        ro.r(
            f'group_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "{group}", label = TRUE, label.size = 5, repel = TRUE)')
        # create a temporary file to save the plot
        with tempfile.NamedTemporaryFile(delete=False, suffix='.png') as temp_file:
            temp_filename = temp_file.name
            # use the file name to create the plot
            ro.r(f'png("{temp_filename}", width = 800, height = 800, res = 150)')
            # print the plot to the device
            ro.r('print(group_plot)')
            # finalize the plot
            ro.r('dev.off()')
        # read the plot file and convert to base64
        with open(temp_filename, "rb") as img_file:
            img_data = img_file.read()
            plot_url = base64.b64encode(img_data).decode('utf-8')
        # remove the temporary file after processing
        os.remove(temp_filename)
        # return the plot as a base64-encoded string to embed in HTML
        return f"data:image/png;base64,{plot_url}"
    except Exception as e:
        return None

# define a function to generate the differential expression table based on the annotation
def generate_table(annotation):
    pandas2ri.activate()
    # load the required seurat library
    ro.r('library(Seurat)')
    ro.r('library(SeuratDisk)')
    # load in the seurat obj
    global seurat_obj
    # r code to extract the markers
    ro.r(f'''
    markers <- FindMarkers(seurat_obj, ident.1 = "{annotation}", group.by = "cell.type.9.long")
    markers <- as.data.frame(markers)
    markers$gene <- rownames(markers)  # Add row names as a new column
    rownames(markers) <- NULL
    markers <- markers[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
    ''')
    # convert result into a pandas dataframe
    with (ro.default_converter + pandas2ri.converter).context():
        markers_df = ro.conversion.get_conversion().rpy2py(ro.r('markers'))
        markers_df.columns = markers_df.columns.str.strip()
        markers_df.columns = ['gene', 'avg_log2FC', 'pct1', 'pct2', 'p_val_adj']
        print("Markers DataFrame:", markers_df)
    return markers_df

# define a function to generate the differential expression table based on the subcluster
def generate_subcluster_table(subcluster):
    pandas2ri.activate()
    # load the required seurat library
    ro.r('library(Seurat)')
    ro.r('library(SeuratDisk)')
    # r code to extract the markers
    ro.r(f'''
    submarkers <- FindMarkers(subcluster_obj, ident.1 = "{subcluster}", group.by = final_celltype)
    submarkers <- as.data.frame(submarkers)
    submarkers$gene <- rownames(submarkers)  # Add row names as a new column
    rownames(submarkers) <- NULL
    submarkers <- submarkers[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
    ''')
    # convert result into a pandas dataframe
    with (ro.default_converter + pandas2ri.converter).context():
        markers_df = ro.conversion.get_conversion().rpy2py(ro.r('submarkers'))
        markers_df.columns = markers_df.columns.str.strip()
        markers_df.columns = ['gene', 'avg_log2FC', 'pct1', 'pct2', 'p_val_adj']
        print("Markers DataFrame:", markers_df)
    return markers_df

# define a function to generate the plot for timepoint
def timepoint_plot():
    pandas2ri.activate()
    # load the seurat library
    ro.r('library(Seurat)')
    ro.r('library(SeuratDisk)')
    ro.r('library(ggplot2)')
    # subset the object based on the selected cell_type
    ro.r(f'timepoint_plot <- DimPlot(subcluster_obj, reduction = "umap", group.by = "timepoint", repel = TRUE)')
    # create a temporary file to save the plot
    with tempfile.NamedTemporaryFile(delete=False, suffix='.png') as temp_file:
        temp_filename = temp_file.name
        # use the file name to create the plot
        ro.r(f'png("{temp_filename}", width = 800, height = 800, res = 150)')
        # print the plot to the device
        ro.r('print(timepoint_plot)')
        # finalize the plot
        ro.r('dev.off()')
    # read the plot file and convert to base64
    with open(temp_filename, "rb") as img_file:
        img_data = img_file.read()
        plot_url = base64.b64encode(img_data).decode('utf-8')
    # remove the temporary file after processing
    os.remove(temp_filename)
    # return the plot as a base64-encoded string to embed in HTML
    return f"data:image/png;base64,{plot_url}"

# define a function to generate the plot for timepoint
def timepoint_plot_violin():
    pandas2ri.activate()
    # load the r libraries
    ro.r('library(Seurat)')
    ro.r('library(SeuratDisk)')
    ro.r('library(ggplot2)')
    # get the gene from the session 
    gene = session.get('subcluster_gene')
    if gene == None:
        return None
    # subset the object based on the selected cell_type
    ro.r(f'timepoint_plot <- VlnPlot(subcluster_obj, features = "{gene}", group.by = "timepoint")')
    # create a temporary file to save the plot
    with tempfile.NamedTemporaryFile(delete=False, suffix='.png') as temp_file:
        temp_filename = temp_file.name
        # use the file name to create the plot
        ro.r(f'png("{temp_filename}", width = 800, height = 800, res = 150)')
        # print the plot to the device
        ro.r('print(timepoint_plot)')
        # finalize the plot
        ro.r('dev.off()')
    # read the plot file and convert to base64
    with open(temp_filename, "rb") as img_file:
        img_data = img_file.read()
        plot_url = base64.b64encode(img_data).decode('utf-8')
    # remove the temporary file after processing
    os.remove(temp_filename)
    # return the plot as a base64-encoded string to embed in HTML
    return f"data:image/png;base64,{plot_url}"

# define a function to generate the plot for a specific cell type
def sub_cluster_plot(cell_type):
    pandas2ri.activate()
    # load in the seurat libraries
    ro.r('library(Seurat)')
    ro.r('library(SeuratDisk)')
    ro.r('library(ggplot2)')
    # clear any subcluster_obj from the environment to prevent conflicts
    ro.r('if (exists("subcluster_obj")) rm(subcluster_obj, envir = .GlobalEnv)')
    # read in the global variables
    global Hepatocyte
    global Biliary_Epithelial_Cell
    global Endothelial_Cell
    global Fibroblast
    global Macrophage
    global Monocyte
    global Erythrocyte
    global Neuron
    # format the cell type correctly
    cell_type = cell_type.strip()  
    # get the global variable name corresponding to the subtype
    global_cell_type_name = re.sub(r"\s+", "_", cell_type)
    # check if the global variable already exists as a seurat obj (has already been loaded)
    if global_cell_type_name in globals() and globals()[global_cell_type_name] is not None:
        # if yes, get the object
        subcluster_obj = globals()[global_cell_type_name]
        # assign as a variable in R
        ro.globalenv['subcluster_obj'] = subcluster_obj
    else:
        # get the file name for the selected cell type
        data_folder = "subclusterData/"
        temp_name = cell_type.replace(" ", "")  # ensure there are no leading/trailing spaces
        file_name = f"{temp_name}SubCluster.h5Seurat"
        # construct the full file path
        file_path = os.path.abspath(os.path.join(data_folder, file_name))
        # load in the seurat obj
        ro.r(f'subcluster_obj <- LoadH5Seurat("{file_path}")')
        ro.r('assign("subcluster_obj", subcluster_obj, envir = .GlobalEnv)')
        # load the seurat obj into python
        subcluster_obj = ro.r('subcluster_obj')
        # store as a global variable
        globals()[global_cell_type_name] = subcluster_obj
    # get the cell type in the right format to retrieve the appropriate group by var
    ro.r(f'temp_celltype <- gsub(" ", "_", "{cell_type}")')
    ro.r(f'final_celltype <- paste0(temp_celltype, "_subcluster")')
    # store the subcluster name into a python session obj
    final_cell_type = ro.r(f'final_celltype')
    session['final_cell_type'] = final_cell_type[0]
    # create the plot
    ro.r(f'cluster_plot <- DimPlot(subcluster_obj, reduction = "umap", group.by = final_celltype, label = TRUE)')
    # create a temporary file to save the plot
    with tempfile.NamedTemporaryFile(delete=False, suffix='.png') as temp_file:
        temp_filename = temp_file.name
        # use the file name to create the plot
        ro.r(f'png("{temp_filename}", width = 800, height = 800, res = 150)')
        # print the plot to the device
        ro.r('print(cluster_plot)')
        # finalize the plot
        ro.r('dev.off()')
    # read the plot file and convert to base64
    with open(temp_filename, "rb") as img_file:
        img_data = img_file.read()
        plot_url = base64.b64encode(img_data).decode('utf-8')
    # remove the temporary file after processing
    os.remove(temp_filename)
    # return the plot as a base64-encoded string to embed in HTML
    return f"data:image/png;base64,{plot_url}"
 
# create a function to make the gene plot
def gene_plot_subcluster(gene):
    pandas2ri.activate()
    # import the r libraries
    ro.r('library(Seurat)')
    ro.r('library(SeuratDisk)')
    ro.r('library(ggplot2)')
    # store the chosen gene in the session
    session['subcluster_gene'] = gene
    # check if the selected gene exists
    gene_check = ro.r(f'"{gene}" %in% rownames(subcluster_obj)')
    if not gene_check[0]:
        return None # return nothing if no
    # create the plot
    ro.r(f'gene_plot <- FeaturePlot(subcluster_obj, features = "{gene}", order=TRUE, min.cutoff = "q10", repel = TRUE)')
    # create a temporary file to save the plot
    with tempfile.NamedTemporaryFile(delete=False, suffix='.png') as temp_file:
        temp_filename = temp_file.name
        # use the file name to create the plot
        ro.r(f'png("{temp_filename}", width = 800, height = 800, res = 150)')
        # print the plot to the device
        ro.r('print(gene_plot)')
        # finalize the plot
        ro.r('dev.off()')
        # read the plot file and convert to base64
    with open(temp_filename, "rb") as img_file:
        img_data = img_file.read()
        plot_url = base64.b64encode(img_data).decode('utf-8')
        # remove the temporary file after processing
    os.remove(temp_filename)
    # return the plot as a base64-encoded string to embed in HTML
    return f"data:image/png;base64,{plot_url}"

# create the login route
@app.route('/login', methods=['GET', 'POST'])
def login():
    # set the error to none
    error = None
    # initialize the login form
    form = LoginForm()
    # check validation
    if form.validate_on_submit():
        # extract the usernames and see if they equal the env variables
        if form.username.data != USERNAME or form.password.data != PASSWORD:
            error = 'Invalid Credentials. Please try again.'
        else:
            session['logged_in'] = True
            return redirect(url_for('home'))
    return render_template('login.html', form=form, error=error)

# create the index route
@app.route('/', methods=['GET'])
def index():
    return redirect(url_for('login'))

# create the home route
@app.route('/home', methods=['GET', 'POST'])
def home():
    # ensure the user is logged in
    if not session.get('logged_in'):
        return redirect(url_for('login'))
    # define the global seurat obj
    global seurat_obj
    # create form instances
    fileform = UploadFileForm()
    annotationform = AnnotationForm()
    groupform = GroupForm()
    geneform = GeneForm()
    groupsubplotform = GroupSubplotForm()
    genesubplotform = GeneSubplotForm()
    timepointform = TimepointForm()
    subclusterform = SubClusterForm()
    # list existing files in the uploads directory
    existing_files = os.listdir(app.config['UPLOAD_FILES_DEST'])
    # NOTE CURRENT COMMENT OUT UPLOAD CODE - WILL BE CONTROLLED BY VM OWNER
    # new file being uploaded
    if request.method == 'POST':
	'''
        if fileform.validate_on_submit():
            # get the file value from the form
            file = fileform.file.data
            filepath = os.path.join(app.config['UPLOAD_FILES_DEST'], file.filename)
            file.save(filepath)
            flash(f'"{file.filename}" uploaded successfully', 'success')
            # reset annotations and selected group in the session after uploading
            session['filepath'] = filepath
            # load seurat
            ro.r('library(Seurat)')
            ro.r('library(SeuratDisk)')
            ro.r(f"seurat_obj <- LoadH5Seurat('{filepath}')")
            # store in the global variable
            seurat_obj = ro.r('seurat_obj')
            # set empty list to store col names
            metadata_columns = []
            annotcols = []
            celltypecols = []
            # get the col names for the umap plot
            with localconverter(pandas2ri.converter):
                meta_data_cols = ro.r(f'colnames(seurat_obj@meta.data)')
                metadata_columns = [str(col) for col in meta_data_cols]
            # get the annotation names for the table
            with localconverter(pandas2ri.converter):
                unique_values = ro.r(f'levels(seurat_obj@meta.data[["cell.type.9.long"]])')
                annotcols = [str(val) for val in unique_values]
            # get the cols for the cell types
            with localconverter(pandas2ri.converter):
                unique_values = ro.r(f'levels(seurat_obj@meta.data[["cell.type.9.long"]])')
                celltypecols = [str(val) for val in unique_values]
            # save the values in the session
            session['ann_cols'] = annotcols
            session['obs_columns'] = metadata_columns
            session['celltype_cols'] = celltypecols
            return redirect(url_for('home'))  # redirect to avoid resubmission
	    '''
        # handle case where wanted file already exists
        if 'existing_file' in request.form:
            # retrieve the selected file
            selected_file = request.form['existing_file']
            if selected_file:
                session['filepath'] = 'uploads/' + selected_file
                filepath = session['filepath']
                flash(f'File "{selected_file}" loaded successfully.', 'success')
                # load seurat
                ro.r('library(Seurat)')
                ro.r('library(SeuratDisk)')
                ro.r(f"seurat_obj <- LoadH5Seurat('{filepath}')")
                # store in the global variable
                seurat_obj = ro.r('seurat_obj')
                # set lists to store the col names
                metadata_columns = []
                annotcols = []
                celltypecols = []
                # get the col names for the umap plot
                with localconverter(pandas2ri.converter):
                    meta_data_cols = ro.r(f'colnames(seurat_obj@meta.data)')
                    metadata_columns = [str(col) for col in meta_data_cols]
                # get the annotation names for the table
                with localconverter(pandas2ri.converter):
                    unique_values = ro.r(f'levels(seurat_obj@meta.data[["cell.type.9.long"]])')
                    annotcols = [str(val) for val in unique_values]
                # get the cols for the cell types
                with localconverter(pandas2ri.converter):
                    unique_values = ro.r(f'levels(seurat_obj@meta.data[["cell.type.9.long"]])')
                    celltypecols = [str(val) for val in unique_values]
                # save the cols to the sessions
                session['ann_cols'] = annotcols
                session['obs_columns'] = metadata_columns
                session['celltype_cols'] = celltypecols
                return redirect(url_for('home'))

    # populate group and annotation forms based on session-stored names
    if 'obs_columns' in session:
        groupform.group.choices = [(col, col) for col in session['obs_columns']]
    if 'ann_cols' in session:
        annotationform.annotation.choices = [(annotation, annotation) for annotation in session['ann_cols']]
    # populate the cell type forms as well
    if 'celltype_cols' in session:
        groupsubplotform.group.choices = [(col, col) for col in session['celltype_cols']]
    # populate the subcluster form
    if 'subcluster_cols' in session:
        subclusterform.cluster.choices = [(col, col) for col in session['subcluster_cols']]


    return render_template('index.html', fileform=fileform, annotationform=annotationform, groupform=groupform,
                           geneform=geneform, existing_files=existing_files, groupsubplotform=groupsubplotform, genesubplotform=genesubplotform,
                           timepointform=timepointform, subclusterform=subclusterform)

# create the gene plot (overall, not subcluster)
@app.route('/plot_gene', methods=['POST'])
def plot_gene():
    # create the form instances
    groupform = GroupForm()
    geneform = GeneForm()
    #fileform = UploadFileForm()
    annotationform = AnnotationForm()

    # populate the form columns to maintain across submissions
    if 'obs_columns' in session:
        groupform.group.choices = [(col, col) for col in session['obs_columns']]
    if 'ann_cols' in session:
        annotationform.annotation.choices = [(annotation, annotation) for annotation in session['ann_cols']]
    # populate the selected choices to maintain across submissions
    if 'selected_group' in session:
        groupform.group.data = session['selected_group']
    if 'annotation' in session:
        annotationform.annotation.data = session['annotation']

    # validate the form submission
    if geneform.validate_on_submit():
        # get the gene of interest and create the plot
        gene_of_interest = geneform.gene.data
        session['gene_of_interest'] = gene_of_interest
        gene_plot = plot_gene_umap(gene_of_interest)

        if gene_plot is None:
            #if gene is not found, flash an error message
            return jsonify({"error": f"Gene '{gene_of_interest}' not found in the dataset."})

        # return JSON response
        return jsonify(gene_plot=gene_plot)

    return jsonify(error="Invalid form submission"), 400  #return error if form is invalid

# create the route to plot the group plot
@app.route('/plot_group', methods=['POST'])
def plot_group():
    # define the form instances
    groupform = GroupForm()
    geneform = GeneForm()
    #fileform = UploadFileForm()
    annotationform = AnnotationForm()

    # populate the form choices to maintain across submissions
    if 'obs_columns' in session:
        groupform.group.choices = [(col, col) for col in session['obs_columns']]
    if 'ann_cols' in session:
        annotationform.annotation.choices = [(annotation, annotation) for annotation in session['ann_cols']]

    # ensure the gene plot is maintained across submissions
    gene_plot = None
    if 'gene_plot' in request.form:
        gene_plot = request.form['gene_plot']

    # ensure the users choices are maintained across submissions
    if 'gene_of_interest' in session:
        geneform.gene.data = session['gene_of_interest']
    if 'annotation' in session:
        annotationform.annotation.data = session['annotation']

    # validate the form submission
    if groupform.validate_on_submit():
        # create the group plot with the users selected group
        selected_group = groupform.group.data
        session['selected_group'] = selected_group
        group_plot = plot_group_umap(selected_group)

        # if option selected doesnt work for the umap
        if group_plot is None:
            return jsonify({"error": f"Not an informative UMAP option. Please select a different group."})

        # Return JSON response
        return jsonify(group_plot=group_plot)

    return jsonify(error="Invalid form submission"), 400  # Return error if form is invalid

# define the route to make the DE table
@app.route('/make_table', methods=['POST'])
def make_table():
    # create the annotation form instance
    annotationform = AnnotationForm()
    # populate the form choices
    if 'ann_cols' in session:
        annotationform.annotation.choices = [(annotation, annotation) for annotation in session['ann_cols']]
    # validate the form submission
    if annotationform.validate_on_submit():
        # submit the annotation and generate the table
        annotation = annotationform.annotation.data
        session['annotation'] = annotation
        # instead of redirecting return a JSON response
        return jsonify({'status': 'success', 'message': 'Annotation saved'})
    else:
        print("Form validation failed:", annotationform.errors)
        return jsonify({'status': 'error', 'message': 'Form validation failed'}), 400

# define a route to create the subcluster table
@app.route('/make_subcluster_table', methods=['POST'])
def make_subcluster_table():
    # create the form instance
    subclusterform = SubClusterForm()
    # populate the forms choices
    if 'subcluster_cols' in session:
        subclusterform.cluster.choices = [(col, col) for col in session['subcluster_cols']]
    # validate the form submission
    if subclusterform.validate_on_submit():
        # create the table with the users choice
        subcluster = subclusterform.cluster.data
        session['subcluster'] = subcluster
        # instead of redirecting return a JSON response
        return jsonify({'status': 'success', 'message': 'Annotation saved'})
    else:
        print("Form validation failed:", subclusterform.errors)
        return jsonify({'status': 'error', 'message': 'Form validation failed'}), 400

# create the route to generate the table data
@app.route('/api/table_data', methods=['GET'])
def table_data():
    # ensure the user is logged in
    if not session.get('logged_in'):
        return jsonify([])

    # get the selected annotation from the session
    annotation = session.get('annotation')
    filepath = session.get('filepath')

    if annotation and filepath:
        # generate table data based on the annotation
        table_data = generate_table(annotation)
        # convert df to list of dictionaries for JSON response
        json_data = table_data.to_dict(orient='records')
        # save the current table to the session
        return jsonify(json_data)
    else:
        return jsonify([])  # return an empty JSON response if no data

# define the route to generate the subcluster table data
@app.route('/api/subcluster_table_data', methods=['GET'])
def subcluster_table_data():
    # ensure the user is logged in
    if not session.get('logged_in'):
        return jsonify([])

    # get the selected cluster from the session
    subcluster = session.get('subcluster')

    if subcluster:
        # generate table data based on the annotation
        table_data = generate_subcluster_table(subcluster)
        # convert df to list of dictionaries for JSON response
        json_data = table_data.to_dict(orient='records')
        # save the current table to the session
        return jsonify(json_data)
    else:
        return jsonify([])  # return an empty JSON response if no data

# define a route for the group plot
@app.route('/plot_group_subcluster', methods=['POST'])
def plot_group_subcluster():
    # create the form instances
    groupsubplotform = GroupSubplotForm()
    subclusterform = SubClusterForm()

    # populate the form choices
    if 'celltype_cols' in session:
        groupsubplotform.group.choices = [(col, col) for col in session['celltype_cols']]

    # validate the form submission
    if groupsubplotform.validate_on_submit():
        # generate the plot with the selected cell type
        selected_group = groupsubplotform.group.data
        session['selected_celltype'] = selected_group 
        cluster_cols = []
        group_subcluster_plot = sub_cluster_plot(selected_group)
        cell_type = session.get('final_cell_type')

        # if the seurat obj already exists to aviod reloading
        if ro.r('exists("subcluster_obj", envir = .GlobalEnv)'):
            # populate the form for the DE table choices with the cell type clusters
            with localconverter(pandas2ri.converter):
                unique_values = ro.r(f'levels(subcluster_obj@meta.data[["{cell_type}"]])')
                print(unique_values)
                cluster_cols = [str(val) for val in unique_values]
                session['subcluster_cols'] = cluster_cols  # Save the columns in the session
                print(f"Subcluster columns: {session.get('subcluster_cols')}")

        return jsonify(group_subcluster_plot=group_subcluster_plot, cluster_cols=cluster_cols)
    return jsonify(error="Invalid form submission"), 400

# create the route to plot the gene subcluster plot
@app.route('/plot_gene_subcluster', methods=['POST'])
def plot_gene_subcluster():
    # create the form instance
    genesubplotform = GeneSubplotForm()

    # validate the form submission and create the plot with selected gene
    if genesubplotform.validate_on_submit():
        gene_of_interest = genesubplotform.gene.data
        session['gene_of_interest_subcluster'] = gene_of_interest
        gene_subcluster_plot = gene_plot_subcluster(gene_of_interest)
        return jsonify(gene_subcluster_plot=gene_subcluster_plot)
    return jsonify(error="Invalid form submission"), 400

# create the route to plot the timepoint subcluster plot
@app.route('/plot_timepoint_subcluster', methods=['POST'])
def plot_timepoint_subcluster():
    # create the form instance
    timepointform = TimepointForm()

    # validate the form on submission
    if timepointform.validate_on_submit():
        # plot either the umap or violin plot depending on user selection
        response = timepointform.timepoint.data
        if response == 'umap' or response =='Umap':
            t_plot = timepoint_plot()
            return jsonify(t_plot=t_plot)
        else:
            t_plot = timepoint_plot_violin()
            if t_plot is None:
                return jsonify({"error": f"Please select a gene first"})
            return jsonify(t_plot=t_plot)
    return jsonify(error="Invalid form submission"), 400

# define the logout route
@app.route('/logout')
def logout():
    # clear the session and global variables
    session.clear()
    global seurat_obj
    global Hepatocyte
    global Biliary_Epithelial_Cell
    global Endothelial_Cell
    global Fibroblast
    global Macrophage
    global Monocyte
    global Erythrocyte
    global Neuron
    seurat_obj = None
    Hepatocyte = None
    Biliary_Epithelial_Cell = None
    Endothelial_Cell = None
    Hepatic_Stellate_Cell = None
    Fibroblast = None
    Macrophage = None
    Monocyte = None
    Erythrocyte = None
    Neuron = None
    return redirect(url_for('login'))


if __name__ == '__main__':
    app.run(debug=True, threaded=False, use_reloader=False)
