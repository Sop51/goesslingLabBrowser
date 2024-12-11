// Spinner function for upload
function showSpinner() {
    const spinner = document.getElementById('spinner');
    const uploadButton = document.getElementById('uploadButton');

    // Show the spinner and hide the button
    spinner.style.display = 'inline-block';
    uploadButton.style.display = 'none';
}

// Rendering the subcluster annotation table
$(document).ready(function() {
    // Initially hide the table until the necessary data is ready
    $('#subcluster-table').hide();

    // Function to initialize the DataTable
    function initializeTable(data) {
        $('#subcluster-table').show();
        // Destroy any existing DataTable instance before initializing a new one
        if ($.fn.dataTable.isDataTable('#subcluster-table')) {
            $('#subcluster-table').DataTable().clear().destroy();
        }
        // Initialize the DataTable after data is fetched
        var dataTable = $('#subcluster-table').DataTable({
            "data": data,  // Load the data directly into the table
            "scrollY": "75vh",
            "scrollCollapse": true,
            "paging": true,
            "searching": true,
            "lengthMenu": [10, 25, 50, 100],
            "pageLength": 50,
            "columns": [
                { "data": "gene" },
                { "data": "avg_log2FC" },
                { "data": "pct1" },
                { "data": "pct2" },
                { "data": "p_val_adj" }
            ],
            "autoWidth": false, // Disable auto column width
            "columnDefs": [
                { "width": "20%", "targets": [0, 1, 2, 3, 4] }
            ]
        });

        // Show the table now that it is ready
        dataTable.draw();
    }

    // Fetch table data after the form is submitted
    function fetchTableData() {
        $.ajax({
            url: $('#subcluster-table').data('ajax-url'),  // Data URL from the table
            type: 'GET',
            success: function(data) {
                console.log("Table data received:", data);  // Ensure the data is logged correctly
                if (Array.isArray(data)) {
                    // Initialize table after fetching data and making it visible
                    setTimeout(function() {
                        initializeTable(data);  // Initialize the table after a slight delay
                    }, 50);  // Small delay to ensure table is visible
                } else {
                    console.error("Invalid data format:", data);  // More detailed error
                }
            },
            error: function(xhr) {
                console.error('Error fetching table data:', xhr.responseText);
                alert('Error fetching table data: ' + xhr.responseText);
            }
        });
    }

    // Handle form submission
    $('#cluster-form').on('submit', function(e) {
        e.preventDefault();  // Prevent default form submission
        var $submitButton = $(this).find('button[type="submit"]');
        $submitButton.prop('disabled', true);  // Disable submit button to prevent multiple submits

        $.ajax({
            url: $('#cluster-form').data('ajax-url'),  // Use the form's action URL
            type: 'POST',  // Assuming the form is submitting data via POST
            data: $(this).serialize(),
            success: function(response) {
                console.log("Form submitted successfully. Response:", response);
                // After the form submission, fetch the table data
                fetchTableData();
                $submitButton.prop('disabled', false);  // Enable submit button again
            },
            error: function(xhr) {
                alert('Error: ' + xhr.responseText);
                $submitButton.prop('disabled', false);  // Enable submit button in case of error
            }
        });
    });
});

// Rendering the annotation table
$(document).ready(function() {
    // Initially hide the table until the necessary data is ready
    $('#annotation-table').hide();

    // Function to initialize the DataTable
    function initializeTable(data) {
        $('#annotation-table').show();
        // Initialize the DataTable after data is fetched
        var dataTable = $('#annotation-table').DataTable({
            "data": data,  // Load the data directly into the table
            "scrollY": "75vh",
            "scrollCollapse": true,
            "paging": true,
            "searching": true,
            "lengthMenu": [10, 25, 50, 100],
            "pageLength": 50,
            "columns": [
                { "data": "gene" },
                { "data": "avg_log2FC" },
                { "data": "pct1" },
                { "data": "pct2" },
                { "data": "p_val_adj" }
            ],
            "autoWidth": false, // Disable auto column width
            "columnDefs": [
                { "width": "20%", "targets": [0, 1, 2, 3, 4] }
            ]
        });

        // Show the table now that it is ready
        dataTable.draw();
    }

    // Fetch table data after the form is submitted
    function fetchTableData() {
        $.ajax({
            url: $('#annotation-table').data('ajax-url'),  // Data URL from the table
            type: 'GET',
            success: function(data) {
                console.log("Table data received:", data);  // Ensure the data is logged correctly
                if (Array.isArray(data)) {
                    // Initialize table after fetching data and making it visible
                    setTimeout(function() {
                        initializeTable(data);  // Initialize the table after a slight delay
                    }, 50);  // Small delay to ensure table is visible
                } else {
                    console.error("Invalid data format:", data);  // More detailed error
                }
            },
            error: function(xhr) {
                console.error('Error fetching table data:', xhr.responseText);
                alert('Error fetching table data: ' + xhr.responseText);
            }
        });
    }

    // Handle form submission
    $('#annotation-form').on('submit', function(e) {
        e.preventDefault();  // Prevent default form submission
        var $submitButton = $(this).find('button[type="submit"]');
        $submitButton.prop('disabled', true);  // Disable submit button to prevent multiple submits

        $.ajax({
            url: $('#annotation-form').data('ajax-url'),  // Use the form's action URL
            type: 'POST',  // Assuming the form is submitting data via POST
            data: $(this).serialize(),
            success: function(response) {
                console.log("Form submitted successfully. Response:", response);
                // After the form submission, fetch the table data
                fetchTableData();
                $submitButton.prop('disabled', false);  // Enable submit button again
            },
            error: function(xhr) {
                alert('Error: ' + xhr.responseText);
                $submitButton.prop('disabled', false);  // Enable submit button in case of error
            }
        });
    });
});



// Handle gene plot submission
$('#gene-form').on('submit', function(e) {
    e.preventDefault();  // Prevent the default form submission

    // Retrieve the AJAX URL from the data attribute
    var ajaxUrl = $('#gene-plot-url').data('ajax-url');

    $.ajax({
        url: ajaxUrl,  // Use the variable defined above
        type: 'POST',
        data: $(this).serialize(),
        success: function(response) {
            // Check if response contains error
            if (response.error) {
                alert('Error: ' + response.error);
            } else {
                $('#gene-plot-container').html('<img src="' + response.gene_plot + '" alt="Gene UMAP Plot" class="img-fluid">');
            }
        },
        error: function(xhr) {
            if (xhr.status === 404) {
                alert("Gene not found. Please try a different gene.");
            } else {
                alert('An error occurred: ' + xhr.responseText);
            }
        }
    });
});


// Handle group plot submission
$('#group-form').on('submit', function(e) {
    e.preventDefault();  // Prevent the default form submission

    // Retrieve the AJAX URL from the data attribute
    var ajaxUrl = $('#group-plot-url').data('ajax-url');

    $.ajax({
        url: ajaxUrl,  // Use the variable defined above
        type: 'POST',
        data: $(this).serialize(),
        success: function(response) {
            // Check if response contains error
            if (response.error) {
                alert('Error: ' + response.error);
            } else {
                $('#group-plot-container').html('<img src="' + response.group_plot + '" alt="Group UMAP Plot" class="img-fluid">');
            }
        },
        error: function(xhr) {
            alert('Error: ' + xhr.responseText);
        }
    });
});


// handle group subcluster plot submission
$(document).on('submit', '#celltype-subcluster-form', function(e) {
    console.log($('#celltype-plot-subcluster-url').data('ajax-url'));
    e.preventDefault();  // Prevent the default form submission
    console.log("Form submitted!")
    // Retrieve the AJAX URL from the data attribute
    var ajaxUrl = $('#celltype-plot-subcluster-url').data('ajax-url');

    $.ajax({
        url: ajaxUrl,  // Use the variable defined above
        type: 'POST',
        data: $(this).serialize(),
        success: function(response) {
            // Check if response contains error
            if (response.error) {
                alert('Error: ' + response.error);
            } else {
                $('#celltype-plot-subcluster-container').html('<img src="' + response.group_subcluster_plot + '" alt="Group UMAP Plot" class="img-fluid">');
            }
        },
        error: function(xhr) {
            alert('Error: ' + xhr.responseText);
        }
    });
});

$(document).on('submit', '#gene-subcluster-form', function(e) {
    console.log($('#gene-plot-subcluster-url').data('ajax-url'));
    e.preventDefault();
    console.log("Form submitted!");
    
    var ajaxUrl = $('#gene-plot-subcluster-url').data('ajax-url');
    $.ajax({
        url: ajaxUrl,
        type: 'POST',
        data: $(this).serialize(),
        success: function(response) {
            if (response.error) {
                alert('Error: ' + response.error);
            } else {
                console.log(response.group_subcluster_plot);
                $('#gene-plot-subcluster-container').html('<img src="' + response.gene_subcluster_plot + '" alt="Gene UMAP Plot" class="img-fluid">');
            }
        },
        error: function(xhr) {
            alert('Error: ' + xhr.responseText);
        }
    });
});


// handle timepoint subcluster plot submission
$(document).on('submit', '#timepoint-form', function(e) {
    e.preventDefault();  // Prevent the default form submission

    // Retrieve the AJAX URL from the data attribute
    var ajaxUrl = $('#timepoint-plot-subcluster-url').data('ajax-url');

    $.ajax({
        url: ajaxUrl,  // Use the variable defined above
        type: 'POST',
        data: $(this).serialize(),
        success: function(response) {
            // Check if response contains error
            if (response.error) {
                alert('Error: ' + response.error);
            } else {
                $('#timepoint-plot-container').html('<img src="' + response.t_plot + '" alt="Group UMAP Plot" class="img-fluid">');
            }
        },
        error: function(xhr) {
            alert('Error: ' + xhr.responseText);
        }
    });
});

window.addEventListener('beforeunload', function(event) {
    // Send a POST request to logout the user when the tab or page is closed
    navigator.sendBeacon('/logout', JSON.stringify({logout: true}));
});
