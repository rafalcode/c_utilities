/*
write me a gpx parser in c please.
ChatGPT
Certainly! Below is a basic example of a GPX (GPS Exchange Format) parser in C. This code assumes a simple GPX structure and may need to be adapted to handle more complex GPX files with additional elements and attributes.
*/

#include <stdio.h>
#include <string.h>

#define MAX_BUFFER_SIZE 1024

typedef struct {
    double latitude;
    double longitude;
    double elevation;
} Waypoint;

void parse_gpx(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    char buffer[MAX_BUFFER_SIZE];
    Waypoint waypoint;
    int in_trkpt = 0; // Flag to indicate if inside <trkpt> element

    while (fgets(buffer, MAX_BUFFER_SIZE, file)) {
        if (strstr(buffer, "<trkpt") != NULL) {
            sscanf(buffer, "  <trkpt lat=\"%lf\" lon=\"%lf\">",
                   &waypoint.latitude, &waypoint.longitude);
            in_trkpt = 1;
        } else if (strstr(buffer, "</trkpt>") != NULL) {
            in_trkpt = 0;
            printf("Latitude: %lf, Longitude: %lf\n",
                   waypoint.latitude, waypoint.longitude);
        } else if (in_trkpt && strstr(buffer, "<ele>") != NULL) {
            sscanf(buffer, "    <ele>%lf</ele>", &waypoint.elevation);
        }
    }

    fclose(file);
}

int main() {
    const char *gpx_filename = "sample.gpx";
    parse_gpx(gpx_filename);
    return 0;
}
// This simple parser focuses on extracting latitude, longitude, and elevation information from GPX waypoints. Keep in mind ... etc.
