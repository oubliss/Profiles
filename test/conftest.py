"""Suite-wide pytest configuration."""
import matplotlib

# Headless backend for every test, selected once here instead of inside
# individual tests (matplotlib.use is process-wide, so doing it in one test
# silently changed the backend for whichever tests ran after it).
matplotlib.use('Agg')
