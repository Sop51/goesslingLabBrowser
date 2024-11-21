function showSpinner() {
    const spinner = document.getElementById('spinner');
    const uploadButton = document.getElementById('uploadButton');

    // Show the spinner and hide the button
    spinner.style.display = 'inline-block';
    uploadButton.style.display = 'none';
}
$(document).ready(function() {
    var ajaxUrl = $('#annotation-table').data('ajax-url');
    var dataTable = $('#annotation-table').DataTable({
        "ajax": {
            "url": ajaxUrl,
            "dataSrc": ""
        },
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

    // Handle form submission
    $('#annotation-form').on('submit', function(e) {
        e.preventDefault();  // Prevent default form submission
        console.log("Form submitted");  // Log form submission

        $.ajax({
            url: $('#annotation-form').data('ajax-url'),  // Use the form's action URL
            type: 'POST',  // Assuming the form is submitting data via POST
            data: $(this).serialize(),
            success: function(response) {
                console.log("Form submitted successfully. Response:", response);

                // After the form submission, fetch the table data
                fetchTableData();
            },
            error: function(xhr) {
                alert('Error: ' + xhr.responseText);
            }
        });
    });

    // Function to fetch table data
    function fetchTableData() {
        $.ajax({
            url: $('#annotation-table').data('ajax-url'),  // Data URL from the table
            type: 'GET',
            success: function(data) {
                console.log("Table data received:", data);  // Ensure the data is logged correctly
                if (Array.isArray(data)) {
                    dataTable.clear().rows.add(data).draw();
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
            alert('Error: ' + xhr.responseText);
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


