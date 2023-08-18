// Chatgpt code.
// that last chagpt left out dattime, here with KML I try and get it to parse that oo.
#include <stdio.h>
#include <string.h>

#define MAX_BUFFER_SIZE 1024

typedef struct {
    double latitude;
    double longitude;
    double altitude;
    char datetime[20];
} Placemark;

void parse_kml(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    char buffer[MAX_BUFFER_SIZE];
    Placemark placemark;
    int in_placemark = 0;
    int in_coordinates = 0;
    int in_datetime = 0;

    while (fgets(buffer, MAX_BUFFER_SIZE, file)) {
        if (strstr(buffer, "<Placemark>") != NULL) {
            in_placemark = 1;
        } else if (strstr(buffer, "</Placemark>") != NULL) {
            in_placemark = 0;
            printf("Latitude: %lf, Longitude: %lf, Altitude: %lf, Datetime: %s\n",
                   placemark.latitude, placemark.longitude, placemark.altitude, placemark.datetime);
        } else if (in_placemark) {
            if (strstr(buffer, "<coordinates>") != NULL) {
                in_coordinates = 1;
            } else if (strstr(buffer, "</coordinates>") != NULL) {
                in_coordinates = 0;
            } else if (in_coordinates) {
                sscanf(buffer, "%lf,%lf,%lf",
                       &placemark.longitude, &placemark.latitude, &placemark.altitude);
            } else if (strstr(buffer, "<when>") != NULL) {
                in_datetime = 1;
            } else if (strstr(buffer, "</when>") != NULL) {
                in_datetime = 0;
            } else if (in_datetime) {
                sscanf(buffer, " %19[^\n]", placemark.datetime);
            }
        }
    }

    fclose(file);
}

int main() {
    const char *kml_filename = "sample.kml";
    parse_kml(kml_filename);
    return 0;
}

