# Check if renv packages need to be restored
if (!requireNamespace("yaml", quietly = TRUE) || !requireNamespace("here", quietly = TRUE)) {
  cat("🔄 First-time setup: Installing packages from renv lockfile...\n")
  cat("   (This may take 10-20 minutes)\n\n")
  renv::restore(prompt = FALSE)
  cat("\n✅ Package installation complete!\n\n")
}