#!/bin/bash

# HeartMAP Hugging Face Deployment Script

echo "🚀 Preparing HeartMAP for Hugging Face Spaces deployment..."

# Set variables
HF_USERNAME="Tumo505" 
SPACE_NAME="HeartMAP" 
REPO_URL="https://huggingface.co/spaces/$HF_USERNAME/$SPACE_NAME"

echo "📁 Creating deployment directory..."
mkdir -p huggingface_deployment
cd huggingface_deployment

# Initialize git repository
echo "🔧 Initializing git repository..."
git init
git remote add origin $REPO_URL

# Copy necessary files
echo "📋 Copying files..."
cp ../huggingface_files/.gitignore .
cp ../huggingface_files/app.py .
cp ../huggingface_files/requirements.txt .
cp ../huggingface_files/README.md .
cp ../huggingface_files/config.yaml .
cp -r ../huggingface_files/src .

# Verify required files
echo " Verifying required files..."
required_files=("app.py" "requirements.txt" "README.md" "config.yaml" ".gitignore")

for file in "${required_files[@]}"; do
    if [[ -f "$file" ]]; then
        echo "  ✓ $file"
    else
        echo "   Missing: $file"
        exit 1
    fi
done

# Check source code
if [[ -d "src" ]]; then
    echo "  ✓ src/ directory"
else
    echo "   Missing: src/ directory"
    exit 1
fi

echo "📦 Files ready for deployment:"
ls -la

echo ""
echo "🎯 Next steps for deployment:"
echo "1. Create a new Space at: https://huggingface.co/new-space"
echo "   - Space name: $SPACE_NAME"
echo "   - SDK: Gradio"
echo "   - Hardware: CPU basic (free)"
echo ""
echo "2. Navigate OUT of this project and clone your space:"
echo "   cd .."
echo "   git clone $REPO_URL"
echo "   cd $SPACE_NAME"
echo ""
echo "3. Copy the prepared files:"
echo "   cp -r ../$(basename $(pwd))/huggingface_deployment/* ."
echo ""
echo "4. Commit and push:"
echo "   git add ."
echo "   git commit -m 'Deploy HeartMAP platform'"
echo "   git push"
echo ""
echo "⚠️  IMPORTANT: Clone the HF space OUTSIDE this project directory!"
echo " Your HeartMAP space will be available at:"
echo "   $REPO_URL"

cd ..
