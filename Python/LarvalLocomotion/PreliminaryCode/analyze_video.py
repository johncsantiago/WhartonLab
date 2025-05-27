import cv2
import csv
import random

# Parameters
video_path = "./Converted/output_video.mp4"
bbox_size = 20  # Width and height of each tracking box
points = []

# Mouse callback to collect click coordinates
def click_event(event, x, y, flags, param):
    if event == cv2.EVENT_LBUTTONDOWN:
        points.append((x, y))
        cv2.circle(frame_copy, (x, y), 3, (0, 255, 0), -1)
        cv2.rectangle(frame_copy,
                      (x - bbox_size // 2, y - bbox_size // 2),
                      (x + bbox_size // 2, y + bbox_size // 2),
                      (0, 255, 0), 1)
        cv2.imshow("Click to Select Objects", frame_copy)

# Load video and first frame
cap = cv2.VideoCapture(video_path)
ret, frame = cap.read()
if not ret:
    print("Failed to read video")
    cap.release()
    exit()

# Let user click objects
frame_copy = frame.copy()
cv2.imshow("Click to Select Objects", frame_copy)
cv2.setMouseCallback("Click to Select Objects", click_event)
print("Click on each object you want to track, then press any key...")

cv2.waitKey(0)
cv2.destroyAllWindows()

# Create bounding boxes from click points
bboxes = [(x - bbox_size // 2, y - bbox_size // 2, bbox_size, bbox_size) for (x, y) in points]

# Create a tracker for each bounding box
trackers = []
colors = []  # Store color for each object

for bbox in bboxes:
    tracker = cv2.TrackerCSRT_create()
    tracker.init(frame, tuple(bbox))
    trackers.append(tracker)

    # Assign a random color for each object
    color = [random.randint(0, 255) for _ in range(3)]  # Random RGB color
    colors.append(color)

# Open CSV to write tracking data
with open('tracked_coordinates_click.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Frame", "ObjectID", "X", "Y", "Width", "Height"])

    frame_count = 0
    while True:
        ret, frame = cap.read()
        if not ret:
            break
        frame_count += 1

        # Add frame number in the top-right corner
        cv2.putText(frame, f"Frame: {frame_count}", (frame.shape[1] - 150, 30),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.8, (255, 255, 255), 2)

        # Update and draw each tracker
        for i, tracker in enumerate(trackers):
            success, bbox = tracker.update(frame)
            if success:
                x, y, w, h = [int(v) for v in bbox]
                # Draw the bounding box with a unique color
                cv2.rectangle(frame, (x, y), (x + w, y + h), colors[i], 2)
                cv2.putText(frame, f"Obj {i+1}", (x, y - 10),
                            cv2.FONT_HERSHEY_SIMPLEX, 0.5, colors[i], 1)
                writer.writerow([frame_count, i + 1, x, y, w, h])
            else:
                cv2.putText(frame, f"Obj {i+1} Lost", (20, 30 + 20 * i),
                            cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 0, 255), 2)

        cv2.imshow("Tracking Multiple Objects", frame)

        if cv2.waitKey(1) & 0xFF == ord('q'):
            break

cap.release()
cv2.destroyAllWindows()
